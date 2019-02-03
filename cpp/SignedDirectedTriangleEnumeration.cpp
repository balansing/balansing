#include "mex.h"
#include "matrix.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>

#define POS +1
#define NEG -1

using namespace std;

class Neighbor{
public:
    int n;
    int sign;
    Neighbor(){}
    Neighbor(int _n, int _sign){
        this->n = _n;
        this->sign = _sign;
    }

    bool operator < (const Neighbor& that) const
    {
        return n < that.n;
    }
};

void printVector(int source, const vector<Neighbor>& neighbors){
    if(neighbors.size() > 0){
        mexPrintf("%d => ", source);
        for(const Neighbor& neighbor : neighbors)
            mexPrintf("(%d, %d), ", neighbor.n, neighbor.sign);

        mexPrintf("\n");
    }
}

void countIntersectBySign(const vector<Neighbor>& uN, \
        const vector<Neighbor>& vN, \
        int& pp, int& pm, int& mp, int& mm){
    pp = pm = mp = mm = 0;

    int uCur = 0;
    int vCur = 0;
    int uD = uN.size();
    int vD = vN.size();

    while(uCur < uD && vCur < vD){
        if(uN[uCur].n < vN[vCur].n){
            uCur++;
        }else if(uN[uCur].n > vN[vCur].n){
            vCur++;
        }else{
            if(uN[uCur].sign == POS && vN[vCur].sign == POS) pp++;
            if(uN[uCur].sign == POS && vN[vCur].sign == NEG) pm++;
            if(uN[uCur].sign == NEG && vN[vCur].sign == POS) mp++;
            if(uN[uCur].sign == NEG && vN[vCur].sign == NEG) mm++;
            uCur++;
            vCur++;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double* X_pr = NULL;
    X_pr = mxGetPr(prhs[0]);
    int X_n = mxGetN(prhs[0]);
    int X_m = mxGetM(prhs[0]);

    int max_node_id = -1;
    for(int i = 0; i < X_m; i++){
        int source = X_pr[0*X_m + i];
        int target = X_pr[1*X_m + i];
        int sign   = X_pr[2*X_m + i];

        max_node_id = std::max(max_node_id, source);
        max_node_id = std::max(max_node_id, target);
    }

    vector<vector<Neighbor>> outNeighbors;
    vector<vector<Neighbor>> inNeighbors;

    for(int i = 0; i <= max_node_id; i++){
        outNeighbors.push_back(vector<Neighbor>());
        inNeighbors.push_back(vector<Neighbor>());
    }

    for(int i = 0; i < X_m; i++){
        int source = X_pr[0*X_m + i];
        int target = X_pr[1*X_m + i];
        int sign   = X_pr[2*X_m + i];

        Neighbor outNeighbor(target, sign);
        Neighbor inNeighbor(source, sign);

        outNeighbors[source].push_back(outNeighbor);
        inNeighbors[target].push_back(inNeighbor);
    }

    for(int i = 1; i <= max_node_id; i++){
        std::sort(outNeighbors[i].begin(), outNeighbors[i].end());
        std::sort(inNeighbors[i].begin(), inNeighbors[i].end());

        // printVector(i, outNeighbors[i]);
        // printVector(i, inNeighbors[i]);
    }

    /**
     * type 0 ~ type 7: non-cycle triangle types
     * The type of a non-cycle triangle is defined by the sign of pivot, left, and right edges
     * * pivot edge: the edge whose both end nodes have out-going edges.
     * * left edge: the edge having the same source node with the pivot edge
     * * right edge: the edge whose source node is the target node of the pivot edge
     * type 0: +++
     * type 1: ++-
     * type 2: +-+
     * type 3: +--
     * type 4: -++
     * type 5: -+-
     * type 6: --+
     * type 7: ---
     *
     * type 8 ~ type 11: cycle triangle types
     * The type of a cycle triangle is defined by the number of plus signs.
     * type 8: +++
     * type 9: ++-
     * type 10: +--
     * type 11: ---
     **/
    int numTrianglesByType[12] = {0, };

    // count non-cyclic triangles
    for(int i = 0; i < X_m; i++){
        int source = X_pr[0*X_m + i];
        int target = X_pr[1*X_m + i];
        int sign   = X_pr[2*X_m + i];

        vector<Neighbor> uN = outNeighbors[source];
        vector<Neighbor> vN = outNeighbors[target];
        int pp, pm, mp, mm;
        countIntersectBySign(uN, vN, pp, pm, mp, mm);

        if(sign == POS){
            numTrianglesByType[0] += pp;
            numTrianglesByType[1] += pm;
            numTrianglesByType[2] += mp;
            numTrianglesByType[3] += mm;
        }else{
            numTrianglesByType[4] += pp;
            numTrianglesByType[5] += pm;
            numTrianglesByType[6] += mp;
            numTrianglesByType[7] += mm;
        }
    }

    // count non-cyclic triangles
    for(int i = 0; i < X_m; i++){
        int source = X_pr[0*X_m + i];
        int target = X_pr[1*X_m + i];
        int sign   = X_pr[2*X_m + i];

        vector<Neighbor> uN = inNeighbors[source];
        vector<Neighbor> vN = outNeighbors[target];
        int pp, pm, mp, mm;
        countIntersectBySign(uN, vN, pp, pm, mp, mm);

        if(sign == POS){
            numTrianglesByType[8] += pp;
            numTrianglesByType[9] += pm;
            numTrianglesByType[9] += mp;
            numTrianglesByType[10] += mm;
        }else{
            numTrianglesByType[9] += pp;
            numTrianglesByType[10] += pm;
            numTrianglesByType[10] += mp;
            numTrianglesByType[11] += mm;
        }
    }

    numTrianglesByType[8] /= 3.0;
    numTrianglesByType[9] /= 3.0;
    numTrianglesByType[10] /= 3.0;
    numTrianglesByType[11] /= 3.0;

    int ncols = 12;
    plhs[0] = mxCreateNumericMatrix(1, ncols, mxDOUBLE_CLASS, mxREAL);
    double* counts_pr = mxGetPr(plhs[0]);
    for(int i = 0; i < ncols; i++){
        counts_pr[i] = numTrianglesByType[i];
    }
}
