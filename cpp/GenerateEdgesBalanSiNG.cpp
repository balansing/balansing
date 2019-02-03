#include "mex.h"
#include "matrix.h"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <queue>
#include <random>
#include <algorithm>
#include <cmath>

using namespace std;

std::random_device rd;
std::default_random_engine rng(rd());
std::uniform_real_distribution<double> uni_real(0.0, 1.0);
double* out_p;
double* in_p;

class Result{
public:
    double p;
    double m;
    int u; // source
    int v; // target
};

class Region{
public:
    pair<int, int> S;
    pair<int, int> D;
};

class Edge{
public:
    int s;
    int t;
    int sign;
};

vector<Region> split_quadrant(Region region){
    int S_s = region.S.first;
    int S_t = region.S.second;
    int D_s = region.D.first;
    int D_t = region.D.second;

    int M = (int) floor((S_t - S_s)/2.0 + 1.0);
    Region Q1, Q2, Q3, Q4;
    // Q1 = [ [S_s, S_s + M - 1], [D_s, D_s + M - 1] ];
    Q1.S = pair<int, int>(S_s, S_s + M - 1);
    Q1.D = pair<int, int>(D_s, D_s + M - 1);

    // Q2 = [ [S_s, S_s + M - 1], [D_s + M, D_t] ];
    Q2.S = pair<int, int>(S_s, S_s + M - 1);
    Q2.D = pair<int, int>(D_s + M, D_t);

    // Q3 = [ [S_s + M, S_t], [D_s, D_s + M - 1] ];
    Q3.S = pair<int, int>(S_s + M, S_t);
    Q3.D = pair<int, int>(D_s, D_s + M - 1);

    // Q4 = [ [S_s + M, S_t], [D_s + M, D_t] ];
    Q4.S = pair<int, int>(S_s + M, S_t);
    Q4.D = pair<int, int>(D_s + M, D_t);

    vector<Region> Q;
    Q.push_back(Q1);
    Q.push_back(Q2);
    Q.push_back(Q3);
    Q.push_back(Q4);
    return Q;
}

vector<double> add_noise(vector<double>& a, double u){
    if(u == 0.0) return a;

    vector<double> _a;

    double a11 = a[0];
    double a12 = a[1];
    double a21 = a[2];
    double a22 = a[3];

    double n11 = (-2.0 * u * a11)/(a11 + a22);
    double n12 = u;
    double n21 = u;
    double n22 = (-2.0 * u * a22)/(a11 + a22);

    _a.push_back(a11 + n11);
    _a.push_back(a12 + n12);
    _a.push_back(a21 + n21);
    _a.push_back(a22 + n22);

    return _a;
}

int weighted_random_selection(vector<double>& a){
    //double indices[4] = {0, 1, 2, 3};
    /*double sum_of_weight = 0;
    for(int i = 0; i < a.size(); i++)
        sum_of_weight += a[i];*/

    double p = uni_real(rng);

    for(int i = 0; i < a.size(); i++){
        //mexPrintf("%.4f and %.4f\n", p, a[i]);
        if(p < a[i]) return i;
        p -= a[i];
    }
    return 0;
}

Result edge_generation(int k, vector<vector<double>>& as, \
        vector<vector<double>>& ps, vector<vector<double>>& ms, \
        Region region, double alpha){

    vector<Region> Q = split_quadrant(region);
    int idx = weighted_random_selection(as[k]);
    //mexPrintf("idx: %d\n", idx);

    double s_p = ps[k][idx];
    double s_m = ms[k][idx];
    Region s_Q = Q[idx];

    Result xk;
    if(k == 1){
        xk.p = s_p;
        xk.m = s_m;
        xk.u = s_Q.S.first;
        xk.v = s_Q.D.first;
    }else{
        Result xk_1 = edge_generation(k - 1, as, ps, ms, s_Q, alpha);

        // balance aggregator
        double tmp_p = xk_1.p * s_p + xk_1.m * s_m;
        double tmp_m = xk_1.p * s_m + xk_1.m * s_p;
        //xk.p = tmp_p;
        //xk.m = tmp_m;

        // alpha aggregator
        xk.p = tmp_p + (1.0 - alpha) * abs(tmp_m);
        xk.m = alpha * tmp_m;

        //xk.p = alpha * tmp_p;
        //xk.m = (1.0 - alpha) * tmp_m;

        xk.u = xk_1.u;
        xk.v = xk_1.v;
    }

    return xk;
}

void gamma_split(vector<double>& a, vector<double>& p, vector<double>& m, double r){
    double a11 = a[0];
    double a12 = a[1];
    double a21 = a[2];
    double a22 = a[3];

    p.push_back( (1 - r) * a11 );
    p.push_back( r * a12 );
    p.push_back( r * a21 );
    p.push_back( (1 - r) * a22 );

    m.push_back( -r * a11 );
    m.push_back( -(1 - r) * a12 );
    m.push_back( -(1 - r) * a21 );
    m.push_back( -r * a22 );
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int L = mxGetScalar(prhs[0]);
    double* A_pr = mxGetPr(prhs[1]);
    double* S_pr = mxGetPr(prhs[2]);
    double* D_pr = mxGetPr(prhs[3]);
    double alpha = mxGetScalar(prhs[4]);
    double b = mxGetScalar(prhs[5]);
    int E = (int) mxGetScalar(prhs[6]);
    double gamma = mxGetScalar(prhs[7]);
    int sign_method = mxGetScalar(prhs[8]);
    int seed = mxGetScalar(prhs[9]);
    
    rng.seed(seed);

    int N = mxGetN(prhs[1]);

    vector<double> a;
    for(int i = 0; i < N; i++){
        double ai1 = A_pr[0 * N + i];
        double ai2 = A_pr[1 * N + i];
        a.push_back(ai1);
        a.push_back(ai2);
    }

    vector<double> p, m;
    gamma_split(a, p, m, gamma);

    std::uniform_real_distribution<double> n_uni_real(-b, b);
    std::uniform_real_distribution<double> alpha_uni_real(0, 1);

    vector<vector<double>> as;
    vector<vector<double>> ps;
    vector<vector<double>> ms;
    vector<double> alphas;
    as.push_back(vector<double>());
    ps.push_back(vector<double>());
    ms.push_back(vector<double>());

    for(int i = 1; i <= L; i++){
        double a_noise = n_uni_real(rng);
        
        vector<double> _a = add_noise(a, a_noise);
        vector<double> _p, _m;
        
        gamma_split(_a, _p, _m, gamma);

        as.push_back(_a);
        ps.push_back(_p);
        ms.push_back(_m);

        double _alpha = alpha_uni_real(rng);
        alphas.push_back(_alpha);
    }

    vector<Edge> edges;

    int pos_edges = 0;
    int neg_edges = 0;

    int no_push_count = 0;

    do{
        Region region;
        region.S = pair<int, int>((int) S_pr[0], (int) S_pr[1]);
        region.D = pair<int, int>((int) D_pr[0], (int) D_pr[1]);

        Result ret = edge_generation(L, as, ps, ms, region, alpha);

        double pk = ret.p;
        double mk = abs(ret.m);

        double p = pk;
        double m = mk;

        if(ret.u == ret.v) continue;

        Edge edge;
        edge.s = ret.u;
        edge.t = ret.v;

        double Ppuv = p/(p + abs(m));
        double Pmuv = 1.0 - Ppuv;

        if(sign_method == 1){
            // determinitic sign decision
            if(Ppuv > Pmuv) edge.sign = +1;
            else if(Ppuv < Pmuv) edge.sign = -1;
            else{
                if(uni_real(rng) <= 0.5) edge.sign = +1;
                else edge.sign = -1;
            }
        }else{
            // stochastic sign decision
            if(uni_real(rng) <= Ppuv) edge.sign = +1;
            else edge.sign = -1;
        }
        
    
        edges.push_back(edge);
    }while(edges.size() < E);

    int nrows = edges.size();
    int ncols = 3;
    plhs[0] = mxCreateNumericMatrix(nrows, ncols, mxDOUBLE_CLASS, mxREAL);
    double* newX_pr = mxGetPr(plhs[0]);

    for(int i = 0; i < nrows; i++){
        Edge edge = edges[i];
        newX_pr[i + 0*nrows] = edge.s;
        newX_pr[i + 1*nrows] = edge.t;
        newX_pr[i + 2*nrows] = edge.sign;
    }
}
