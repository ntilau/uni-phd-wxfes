/**************************************************************************
 *
 * function [P,UV,TRI]=nrbSrfRegularEvalIGES(nurbs,umin,umax,nu,vmin,vmax,nv)
 *
 * Evaluates a NURBS surface at parameter values in a regular grid
 * (see e.g. plotNURBS.m)
 *
 * Usage in Matlab:
 *
 * P=nrbSrfRegularEvalIGES(nurbs,umin,umax,nu,vmin,vmax,nv)
 * [P,UV]=nrbSrfRegularEvalIGES(nurbs,umin,umax,nu,vmin,vmax,nv)
 * [P,UV,TRI]=nrbSrfRegularEvalIGES(nurbs,umin,umax,nu,vmin,vmax,nv)
 *
 * Input:
 * nurbs - NURBS structure
 * umin - fist u-value
 * umax - last u-value
 * nu - number of u-values
 * vmin - fist v-value
 * vmax - last v-value
 * nv - number of v-values
 *
 * Output:
 * P - Points (evaluated NURBS)
 * UV - parameter values of P
 * TRI - triangulation of the regular grid
 *
 * c-file can be downloaded for free at
 *
 * http://www.mathworks.com/matlabcentral/fileexchange/13253-iges-toolbox
 *
 * compile in Matlab by using the command  "mex nrbSrfRegularEvalIGES.c"
 *
 * See "help mex" for more information
 *
 * written by Per Bergström 2012-02-15
 * per.bergstrom at ltu.se
 *
 **************************************************************************/

#include <math.h>
#include "mex.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif


/* Input Arguments */

#define	nurbsstructure	prhs[0]
#define umin	prhs[1]
#define	umax	prhs[2]
#define nu	prhs[3]
#define vmin	prhs[4]
#define	vmax	prhs[5]
#define nv	prhs[6]

/* Output Arguments */

#define	evaluated_points	plhs[0]
#define	uvregular	plhs[1]
#define	tri	plhs[2]


/* Sub functions (in folder "mexSourceFiles") */

#include "mexSourceFiles/FindSpan.c"
#include "mexSourceFiles/FindSpanIncSeq.c"
#include "mexSourceFiles/BasisFuns.c"
#include "mexSourceFiles/NURBSsurfaceRegularEval.c"

/* Main function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    int i, j, intu, intv, myint, myint2, myint3;
    int numu, numv;
    double u0, v0, deltau, deltav, param;
    unsigned int *ptri=NULL;
    int dims[2], mtri;
    
    if (nlhs==0 || nrhs<7){
        mexErrMsgTxt("Wrong number of inputs or outputs.");
    }
    
    numu = (int)mxGetScalar(nu);
    numv = (int)mxGetScalar(nv);
    
    if (numu<2 || numv<2){
        mexErrMsgTxt("nu and nv must be >1");
    }

    intu=numu-1;
    intv=numv-1;
    
    u0=mxGetScalar(umin);
    v0=mxGetScalar(vmin);    
    
    deltau=(mxGetScalar(umax)-u0)/intu;
    deltav=(mxGetScalar(vmax)-v0)/intv;
    
    if (deltau<0 || deltav<0){
        mexErrMsgTxt("umin<=umax and vmin<=vmax are not both true");
    }    
    
    evaluated_points = mxCreateDoubleMatrix(3, numu*numv, mxREAL);
    
    NURBSsurfaceRegularEval((int)mxGetPr(mxGetField(nurbsstructure, 0, "order"))[0], (int)mxGetPr(mxGetField(nurbsstructure, 0, "order"))[1], mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), (int)mxGetDimensions(mxGetField(nurbsstructure, 0, "coefs"))[1], (int)mxGetDimensions(mxGetField(nurbsstructure, 0, "coefs"))[2], mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0)), mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1)), u0, deltau, numu, v0, deltav, numv, mxGetPr(evaluated_points));
    
    if(nlhs>1){
        
        uvregular = mxCreateDoubleMatrix(2, numu*numv, mxREAL);
        
        param=u0;
        for (i = 0; i < 2*numu; i+=2){
            mxGetPr(uvregular)[i]=param;
            mxGetPr(uvregular)[i+1]=v0;
            param+=deltau;
        }
        param=v0;
        for (j = 2*numu; j < 2*numv*numu; j+=2*numu){
            param+=deltav;
            mxGetPr(uvregular)[j]=u0;
            mxGetPr(uvregular)[j+1]=param;
        }
        for (j = 2*numu; j < 2*numv*numu; j+=2*numu){
            for (i = 2; i < 2*numu; i+=2){
                mxGetPr(uvregular)[i+j]=mxGetPr(uvregular)[i];
                mxGetPr(uvregular)[i+j+1]=mxGetPr(uvregular)[j+1];
            }
        }
        
        if(nlhs>2){
            
            mtri=2*intu*intv;
            
            dims[0]=mtri;
            dims[1]=3;
            
            tri = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
            ptri = mxGetData(tri);
            
            for (j = 0; j < intv; j++){
                myint=j*intu;
                myint2=j*numu;
                for (i = 0; i < intu; i++){
                    myint3=2*(myint+i);
                            
                    ptri[myint3]=myint2+numu+i+1;
                    ptri[myint3+1]=myint2+i+2;
                    
                    ptri[myint3+mtri]=myint2+i+2;
                    ptri[myint3+1+mtri]=myint2+numu+i+1;
                    
                    ptri[myint3+2*mtri]=myint2+i+1;
                    ptri[myint3+1+2*mtri]=myint2+numu+i+2;
                }
            }
            
        }
    }
    
}
