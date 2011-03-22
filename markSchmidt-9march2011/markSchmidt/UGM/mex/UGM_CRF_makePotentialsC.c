#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Variables */
    int i,n, s, f, e, n1, n2, s1, s2, nInstances,nNodes, nEdges, nNodeFeatures, nEdgeFeatures, *nStates, maxState,
            *nodeMap, *edgeMap, *edgeEnds, sizeNodePot[2], sizeEdgePot[3];
    
    double *w, *Xnode, *Xedge, *nodePot, *edgePot;
    
    
    /* Input */
    w = mxGetPr(prhs[0]);
    Xnode = mxGetPr(prhs[1]);
    Xedge = mxGetPr(prhs[2]);
    nodeMap = mxGetPr(prhs[3]);
    edgeMap = mxGetPr(prhs[4]);
    nStates = mxGetPr(prhs[5]);
    edgeEnds = mxGetPr(prhs[6]);
    i = mxGetScalar(prhs[7]);
    i--;
    
    /* Compute Sizes */
    nNodes = mxGetDimensions(prhs[3])[0];
    nEdges = mxGetDimensions(prhs[6])[0];
    nInstances = mxGetDimensions(prhs[1])[0];
    nNodeFeatures = mxGetDimensions(prhs[1])[1];
    nEdgeFeatures = mxGetDimensions(prhs[2])[1];
    maxState = getMaxState(nStates, nNodes);
    
    /*printf("%d,%d,%d,%d,%d\n", nNodes, nEdges, nNodeFeatures, nEdgeFeatures, maxState);*/
    /*printf("%d,%d\n",nInstances,i);*/
    
    /* Make output */
    sizeNodePot[0] = nNodes;
    sizeNodePot[1] = maxState;
    sizeEdgePot[0] = maxState;
    sizeEdgePot[1] = maxState;
    sizeEdgePot[2] = nEdges;
    plhs[0] = mxCreateNumericArray(2, sizeNodePot, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, sizeEdgePot, mxDOUBLE_CLASS, mxREAL);
    nodePot = mxGetPr(plhs[0]);
    edgePot = mxGetPr(plhs[1]);
    
    
    for(n = 0; n < nNodes; n++) {
        for(s = 0; s < nStates[n]; s++) {
            for(f = 0; f < nNodeFeatures; f++) {
                if(nodeMap[n + nNodes*(s + maxState*f)] > 0) {
                    nodePot[n+nNodes*s] += w[nodeMap[n+nNodes*(s+maxState*f)]-1]*Xnode[i + nInstances*(f + nNodeFeatures*n)];
                }
            }
            nodePot[n+nNodes*s] = exp(nodePot[n+nNodes*s]);
        }
    }
    
    
    for(e = 0; e < nEdges; e++) {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        
        for (s1 = 0; s1 < nStates[n1]; s1++) {
            for (s2 = 0; s2 < nStates[n2]; s2++) {
                for (f = 0; f < nEdgeFeatures; f++) {
                    if (edgeMap[s1 + maxState*(s2 + maxState*(e + nEdges*f))] > 0) {
                        edgePot[s1+maxState*(s2+maxState*e)] += w[edgeMap[s1+maxState*(s2+maxState*(e+nEdges*f))]-1]*Xedge[i + nInstances*(f+nEdgeFeatures*e)];
                    }
                }
                edgePot[s1+maxState*(s2+maxState*e)] = exp(edgePot[s1+maxState*(s2+maxState*e)]);
            }
        }
    }
    
    
}
