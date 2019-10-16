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
#include <cstdlib>

using namespace std;

std::random_device rd;
std::default_random_engine rng(rd());
std::uniform_real_distribution<double> uni_real(0.0, 1.0);


class Edge{
public:
    long s;
    long t;
    int sign;
    Edge(long _s, long _t, int _sign){
        this->s = _s;
        this->t = _t;
        this->sign = _sign;
    }
};


void gamma_split(vector<double>& a, vector<double>& p, vector<double>& m, double r)
{
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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int L = mxGetScalar(prhs[0]);
    double* A_pr = mxGetPr(prhs[1]);
    double* S_pr = mxGetPr(prhs[2]);
    double* D_pr = mxGetPr(prhs[3]);
    double alpha = mxGetScalar(prhs[4]);
    double beta = mxGetScalar(prhs[5]); // noise (gamma)
    int E = (int) mxGetScalar(prhs[6]);
    int seed = (int) mxGetScalar(prhs[7]);
    double gamma = 0.0; // gamma split (depericated)

    rng.seed(seed);


    int N = mxGetN(prhs[1]);

    vector<double> a, p, m;
    for(int i = 0; i < N; i++){
        double ai1 = A_pr[0 * N + i];
        double ai2 = A_pr[1 * N + i];
        a.push_back(ai1);
        a.push_back(ai2);
    }
    gamma_split(a, p, m, gamma);

    std::uniform_real_distribution<double> n_uni_real(-beta, beta);

    vector<vector<double>> as, ps, ms;
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

        //double _alpha = alpha_uni_real(rng);
        //alphas.push_back(_alpha);
    }

    vector<Edge> edges;

    int pos_edges = 0;
    int neg_edges = 0;

    int no_push_count = 0;

    vector<vector<double>> cumuls;
    vector<double> m_ps;
    vector<double> m_ms;
    for(int l = L; l >= 0; l--){
        vector<double> cumul;
        cumuls.push_back(cumul);

        m_ps.push_back(0.0);
        m_ms.push_back(0.0);
    }

    for(int l = L; l > 0; l--){
        cumuls[l].push_back(as[l][0]);
        cumuls[l].push_back(as[l][0] + as[l][1]);
        cumuls[l].push_back(as[l][0] + as[l][1] + as[l][2]);
    }

    do{
        // determine paths from L to 1

        int S_s = 1;
        int S_t = (int) pow(2, L);
        int D_s = 1;
        int D_t = (int) pow(2, L);

        for(int l = L; l >= 1; l--){
            // select region
            double p = uni_real(rng);

            int M = (int) floor((S_t - S_s)/2.0 + 1.0);

            int i = 0;
            if(p < cumuls[l][0]){
                // Q1 (1, 1)
                S_s = S_s; S_t = S_s + M - 1;
                D_s = D_s; D_t = D_s + M - 1;
                i = 0;
            }else if(p < cumuls[l][1]){
                // Q2 (1, 2)
                S_s = S_s; S_t = S_s + M - 1;
                D_s = D_s + M; D_t = D_t;
                i = 1;
            }else if(p < cumuls[l][2]){
                // Q3 (2, 1)
                S_s = S_s + M; S_t = S_t;
                D_s = D_s; D_t = D_s + M - 1;
                i = 2;
            }else{
                // Q4 (2, 2)
                S_s = S_s + M; S_t = S_t;
                D_s = D_s + M; D_t = D_t;
                i = 3;
            }

            m_ps[l] = ps[l][i];
            m_ms[l] = ms[l][i];
        }

        int u = S_s;
        int v = D_t;

        double p_P = m_ps[1];
        double p_M = m_ms[1];

        for(int l = 2; l <= L; l++){
            double tmp_p = p_P * m_ps[l] + p_M * m_ms[l];
            double tmp_m = p_P * m_ms[l] + p_M * m_ps[l];

            p_P = tmp_p + alpha * abs(tmp_m);
            p_M = (1.0 - alpha) * tmp_m;
        }

        //mexPrintf("%d %d\n", u, v);
        double Ppuv = p_P / (p_P + abs(p_M));

        int sign = -1.0;
        if(uni_real(rng) < Ppuv) sign = +1.0;
        edges.push_back(Edge(u, v, sign));
        //
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
