/*=================================================================
 *
 * This is an implementation of the ICP method using a DVG-tree [1].
 * Tukey's biweight-estimator is used to make it more robust and
 * reduce the influence of bad data (Beaton/Tukey) [2].
 *
 * Commpile in Matlab by running: mex -v icpDVG.c
 *
 * Usage:
 *
 * [R,t]=icpDVG(DATA,quality,quasirand,sizerand,randincreasing,maxiter,kTuSq,MODEL,DIRVEC,normal,pdir,sdir,DVGdata)
 *
 * Input:
 *
 * DATA - the data points
 * quality - a vector with quality values corresponding to DATA
 * quasirand - a vector with a pseudo random permumation of the indices of DATA
 * sizerand - size of the random subset of DATA which is used in the iterations
 * randincreasing - the increasing of sizerand in each iteration 
 * maxiter - the number of ICP-iterations
 * kTuSq - the squared value of the parameter in Tukey's biweight-estimator, see [2]
 * MODEL - the model points
 * DIRVEC - an array with linear approximations (output from icpSrfLinRep)
 * normal - normal direction
 * pdir - primary direction
 * sdir - secondary direction
 * DVGdata - the data structure (output from getDVGtree)
 *
 * Output:
 *
 * R - rotation matrix
 * t - translation vector
 *
 * Written by Per Bergström 2012-03-02
 *
 * email: per.bergstrom@ltu.se
 *
 *
 * References
 *
 * [1]
 * Authors: Per Bergström, Ove Edlund, and Inge Söderkvist
 * Title: Repeated surface registration for on-line use
 * Journal: The International Journal of Advanced Manufacturing Technology
 * Cover Date: 2011-05-01
 * Publisher: Springer London
 * Issn: 0268-3768
 * Pages: 677-689
 * Volume: 54
 * Issue: 5
 * Url: http://dx.doi.org/10.1007/s00170-010-2950-6
 * Doi: 10.1007/s00170-010-2950-6
 *
 * [2]
 * Author: Ralf Wolke
 * Title: Iteratively reweighted least squares: A comparison of several single step algorithms for linear models
 * Journal: BIT Numerical Mathematics
 * Cover Date: 1992-09-01
 * Publisher: Springer Netherlands
 * Issn: 0006-3835
 * Pages: 506-524
 * Volume: 32
 * Issue: 3
 * Url: http://dx.doi.org/10.1007/BF02074884
 * Doi: 10.1007/BF02074884
 *
 *=================================================================*/


#include <math.h>
#include "mex.h"


/* Input Arguments */

#define	data	prhs[0]
#define	quality	prhs[1]
#define	randvec	prhs[2]
#define	sizerand	prhs[3]
#define	randincrease	prhs[4]
#define	maxiter	prhs[5]
#define	kTukeySq	prhs[6]
#define	model	prhs[7]
#define	dirvec	prhs[8]
#define	normal	prhs[9]
#define	pdir	prhs[10]
#define	sdir	prhs[11]
#define	DVGdata	prhs[12]


/* Output Arguments */

#define	rotationMatrixOut	plhs[0]
#define	transVecOut     plhs[1]


#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define	caxis	0


double pwr2(double a){
    return a*a;
}

void icpDVG(
        double *dataPoints,
        double *qltyDataPoints,
        unsigned int numDataPoints,
        unsigned int *quasiRndNumbers,
        unsigned int numRndInIters,
        unsigned int rndIncrease,
        unsigned int ICP_iters,
        double sqrdTukeysParam,
        double *modelPoints,
        double *rctApprx,
        double *gridNormalDir,
        double *gridPrimalDir,
        double *gridSecDir,
        mxArray *DVGtree,
        double *rotationMatrix,
        double *translationVector
        ){
    
    unsigned int i, threei, dataIndRndIndex, j, k, gridInd, treeind, ii, indQuasiRndNumbers=0;
    int cloIndTmp, threecloIndTmp, eightcloIndTmp;
    double dataTemp[3], dataOriginalTemp[3], tmpWgh;
    double gridSpace0, *refPnt, gridIndFloat, doubleVar;
    unsigned int n00, n01, n02;
    unsigned int *uintpntr, nT, *Nvertex, *Nvertex_1, *NvertexSq, *Nvertex_1xn, *Nvertex_1xnxn;
    double *dNvertex, *dNvertex_1;
    int *dataStruct0, *dataStruct;
    double datam[3], modelm[3], MM[3], Spx[9], quasum;
    unsigned short int bol;
    double SIGMA[3];
    double SpxTSpx[6];
    double A, B, C;
    double sqrtA23B;
    double x0, f0;
    double SIp, difl;
    double invmat[6];
    double V[9];
    double U[9];
    unsigned int numQuasiRndNumbers_1, n00_1, n01_1, n02_1, n01_1xn00, n02_1xn00xn01, n00xn01;
    double dn00_1, dn01_1, dn02_1;
    
    uintpntr=(unsigned int*)mxGetPr(mxGetCell(DVGtree,0));
    nT=uintpntr[0];
    
    refPnt=mxGetPr(mxGetCell(DVGtree,1));
    
    uintpntr=(unsigned int*)mxGetPr(mxGetCell(DVGtree,2));
    n00=uintpntr[0];
    n01=uintpntr[1];
    n02=uintpntr[2];
    
    gridSpace0=mxGetPr(mxGetCell(DVGtree,3))[0];
    
    Nvertex=(unsigned int*)mxGetPr(mxGetCell(DVGtree,4));
    
    dataStruct0=(int*)mxGetPr(mxGetCell(DVGtree,5));
    
    numQuasiRndNumbers_1=numDataPoints-1;
    
    n00_1=n00-1;
    n01_1=n01-1;
    n02_1=n02-1;
    n01_1xn00=n01_1*n00;
    n02_1xn00xn01=n02_1*n00*n01;
    n00xn01=n00*n01;
    dn00_1=(double)n00_1;
    dn01_1=(double)n01_1;
    dn02_1=(double)n02_1;
    
    Nvertex_1=(unsigned int *)malloc(nT*sizeof(unsigned int));
    NvertexSq=(unsigned int *)malloc(nT*sizeof(unsigned int));
    Nvertex_1xn=(unsigned int *)malloc(nT*sizeof(unsigned int));
    Nvertex_1xnxn=(unsigned int *)malloc(nT*sizeof(unsigned int));
    dNvertex=(double *)malloc(nT*sizeof(double));
    dNvertex_1=(double *)malloc(nT*sizeof(double));
    
    for (i=0;i<nT;i++){
        Nvertex_1[i]=Nvertex[i]-1;
        NvertexSq[i]=Nvertex[i]*Nvertex[i];
        Nvertex_1xn[i]=NvertexSq[i]-Nvertex[i];
        Nvertex_1xnxn[i]=NvertexSq[i]*Nvertex_1[i];
        dNvertex[i]=(double)Nvertex[i];
        dNvertex_1[i]=(double)Nvertex_1[i];
    }
    
    rotationMatrix[0]=1.0;
    rotationMatrix[4]=1.0;
    rotationMatrix[8]=1.0;
    
    if (numRndInIters>=numDataPoints){
        bol=0;
    }
    else {
        bol=1;
        if (numRndInIters<10){
            numRndInIters=10;
        }
    }
    
    
    for (ii=0;ii<ICP_iters;ii++){
        
        Spx[0]=0.0;
        Spx[1]=0.0;
        Spx[2]=0.0;
        Spx[3]=0.0;
        Spx[4]=0.0;
        Spx[5]=0.0;
        Spx[6]=0.0;
        Spx[7]=0.0;
        Spx[8]=0.0;
        
        modelm[0]=0.0;
        modelm[1]=0.0;
        modelm[2]=0.0;
        
        datam[0]=0.0;
        datam[1]=0.0;
        datam[2]=0.0;
        
        quasum=0.0;
        
        for (dataIndRndIndex=0;dataIndRndIndex<numDataPoints;dataIndRndIndex++){
            
            if (bol){
                if (dataIndRndIndex<numRndInIters){
                    i=quasiRndNumbers[indQuasiRndNumbers];
                    if (indQuasiRndNumbers<numQuasiRndNumbers_1){
                        indQuasiRndNumbers++;
                    }
                    else{
                        indQuasiRndNumbers=0;
                    }
                }
                else{
                    break;
                }
            }
            else{
                i=dataIndRndIndex;
            }
            
            tmpWgh=qltyDataPoints[i];
            threei=3*i;
            
            if (tmpWgh>0.0){
                
                dataOriginalTemp[0]=dataPoints[threei];
                dataOriginalTemp[1]=dataPoints[threei+1];
                dataOriginalTemp[2]=dataPoints[threei+2];
                
                dataTemp[0]=rotationMatrix[0]*dataOriginalTemp[0]+rotationMatrix[3]*dataOriginalTemp[1]+rotationMatrix[6]*dataOriginalTemp[2]+translationVector[0];
                dataTemp[1]=rotationMatrix[1]*dataOriginalTemp[0]+rotationMatrix[4]*dataOriginalTemp[1]+rotationMatrix[7]*dataOriginalTemp[2]+translationVector[1];
                dataTemp[2]=rotationMatrix[2]*dataOriginalTemp[0]+rotationMatrix[5]*dataOriginalTemp[1]+rotationMatrix[8]*dataOriginalTemp[2]+translationVector[2];
                
                
                /* Finding the index of the closest model point  */
                
                
                /**** TREE LEVEL 0  ****/
                
                /* normal direction */
                if (caxis){
                    gridIndFloat=(refPnt[2]-dataTemp[2])/gridSpace0;
                }
                else{
                    gridIndFloat=(gridNormalDir[0]*(dataTemp[0]-refPnt[0])+gridNormalDir[1]*(dataTemp[1]-refPnt[1])+gridNormalDir[2]*(dataTemp[2]-refPnt[2]))/gridSpace0;
                }
                
                if (gridIndFloat<1.0){
                    gridInd=0;
                    SIGMA[0]=gridIndFloat;
                }
                else if (gridIndFloat>=n00_1){
                    gridInd=n00_1;
                    SIGMA[0]=gridIndFloat-dn00_1;
                }
                else{
                    gridInd=(unsigned int)gridIndFloat;
                    SIGMA[0]=gridIndFloat-(double)gridInd;
                }
                
                /* primary direction */
                if (caxis){
                    gridIndFloat=(dataTemp[0]-refPnt[0])/gridSpace0;
                }
                else{
                    gridIndFloat=(gridPrimalDir[0]*(dataTemp[0]-refPnt[0])+gridPrimalDir[1]*(dataTemp[1]-refPnt[1])+gridPrimalDir[2]*(dataTemp[2]-refPnt[2]))/gridSpace0;
                }
                
                if (gridIndFloat<1.0){
                    SIGMA[1]=gridIndFloat;
                }
                else if (gridIndFloat>=n01_1){
                    gridInd+=n01_1xn00;
                    SIGMA[1]=gridIndFloat-dn01_1;
                }
                else{
                    k=(unsigned int)gridIndFloat;
                    gridInd+=k*n00;
                    SIGMA[1]=gridIndFloat-(double)k;
                }
                
                /* secondary direction */
                if (caxis){
                    gridIndFloat=(dataTemp[1]-refPnt[1])/gridSpace0;
                }
                else{
                    gridIndFloat=(gridSecDir[0]*(dataTemp[0]-refPnt[0])+gridSecDir[1]*(dataTemp[1]-refPnt[1])+gridSecDir[2]*(dataTemp[2]-refPnt[2]))/gridSpace0;
                }
                
                if (gridIndFloat<1.0){
                    SIGMA[2]=gridIndFloat;
                }
                else if (gridIndFloat>=n02_1){
                    gridInd+=n02_1xn00xn01;
                    SIGMA[2]=gridIndFloat-dn02_1;
                }
                else{
                    k=(unsigned int)gridIndFloat;
                    gridInd+=k*n00xn01;
                    SIGMA[2]=gridIndFloat-(double)k;
                }
                
                treeind=0;
                cloIndTmp=dataStruct0[gridInd];
                
                /**** TREE LEVEL 1,2,...  ****/
                
                while (cloIndTmp<0){
                    
                    SIGMA[0]*=dNvertex[treeind];
                    SIGMA[1]*=dNvertex[treeind];
                    SIGMA[2]*=dNvertex[treeind];
                    
                    /* normal direction */
                    if (SIGMA[0]<1.0){
                        gridInd=0;
                    }
                    else if (SIGMA[0]>=Nvertex_1[treeind]){
                        gridInd=Nvertex_1[treeind];
                        SIGMA[0]-=-dNvertex_1[treeind];
                    }
                    else{
                        gridInd=(unsigned int)SIGMA[0];
                        SIGMA[0]-=(double)gridInd;
                    }
                    
                    /* primary direction */
                    if (SIGMA[1]<1.0){
                        
                    }
                    else if (SIGMA[1]>=Nvertex_1[treeind]){
                        gridInd+=Nvertex_1xn[treeind];
                        SIGMA[1]-=dNvertex_1[treeind];
                    }
                    else{
                        k=(unsigned int)SIGMA[1];
                        gridInd+=k*Nvertex[treeind];
                        SIGMA[1]-=(double)k;
                    }
                    
                    /* secondary direction */
                    if (SIGMA[2]<1.0){
                        
                    }
                    else if (SIGMA[2]>=Nvertex_1[treeind]){
                        gridInd+=Nvertex_1xnxn[treeind];
                        SIGMA[2]-=dNvertex_1[treeind];
                    }
                    else{
                        k=(unsigned int)SIGMA[2];
                        gridInd+=k*NvertexSq[treeind];
                        SIGMA[2]-=(double)k;
                    }
                    
                    dataStruct=(int*)mxGetPr(mxGetCell(DVGtree,6+treeind));
                    cloIndTmp=dataStruct[gridInd-cloIndTmp-1];
                    treeind++;
                    
                }
                
                threecloIndTmp=3*cloIndTmp;
                eightcloIndTmp=8*cloIndTmp;
                
                /* Closest model index is found, MM closest model point */
                
                MM[0]=modelPoints[threecloIndTmp];
                MM[1]=modelPoints[threecloIndTmp+1];
                MM[2]=modelPoints[threecloIndTmp+2];
                
                
                /* Start using the DIRVEC for applying the linear surface approximation */
                
                doubleVar=(dataTemp[0]-MM[0])*rctApprx[eightcloIndTmp+1]+(dataTemp[1]-MM[1])*rctApprx[eightcloIndTmp+2]+(dataTemp[2]-MM[2])*rctApprx[eightcloIndTmp+3];
                if (doubleVar>rctApprx[eightcloIndTmp]){
                    MM[0]+=rctApprx[eightcloIndTmp]*rctApprx[eightcloIndTmp+1];
                    MM[1]+=rctApprx[eightcloIndTmp]*rctApprx[eightcloIndTmp+2];
                    MM[2]+=rctApprx[eightcloIndTmp]*rctApprx[eightcloIndTmp+3];
                }
                else if (doubleVar<(-rctApprx[eightcloIndTmp])){
                    MM[0]-=rctApprx[eightcloIndTmp]*rctApprx[eightcloIndTmp+1];
                    MM[1]-=rctApprx[eightcloIndTmp]*rctApprx[eightcloIndTmp+2];
                    MM[2]-=rctApprx[eightcloIndTmp]*rctApprx[eightcloIndTmp+3];
                }
                else{
                    MM[0]+=doubleVar*rctApprx[eightcloIndTmp+1];
                    MM[1]+=doubleVar*rctApprx[eightcloIndTmp+2];
                    MM[2]+=doubleVar*rctApprx[eightcloIndTmp+3];
                }
                doubleVar=(dataTemp[0]-MM[0])*rctApprx[eightcloIndTmp+5]+(dataTemp[1]-MM[1])*rctApprx[eightcloIndTmp+6]+(dataTemp[2]-MM[2])*rctApprx[eightcloIndTmp+7];
                if (doubleVar>rctApprx[eightcloIndTmp+4]){
                    MM[0]+=rctApprx[eightcloIndTmp+4]*rctApprx[eightcloIndTmp+5];
                    MM[1]+=rctApprx[eightcloIndTmp+4]*rctApprx[eightcloIndTmp+6];
                    MM[2]+=rctApprx[eightcloIndTmp+4]*rctApprx[eightcloIndTmp+7];
                }
                else if (doubleVar<(-rctApprx[eightcloIndTmp+4])){
                    MM[0]-=rctApprx[eightcloIndTmp+4]*rctApprx[eightcloIndTmp+5];
                    MM[1]-=rctApprx[eightcloIndTmp+4]*rctApprx[eightcloIndTmp+6];
                    MM[2]-=rctApprx[eightcloIndTmp+4]*rctApprx[eightcloIndTmp+7];
                }
                else{
                    MM[0]+=doubleVar*rctApprx[eightcloIndTmp+5];
                    MM[1]+=doubleVar*rctApprx[eightcloIndTmp+6];
                    MM[2]+=doubleVar*rctApprx[eightcloIndTmp+7];
                }
                
                doubleVar=pwr2(dataTemp[0]-MM[0])+pwr2(dataTemp[1]-MM[1])+pwr2(dataTemp[2]-MM[2]);
                
                if (doubleVar<sqrdTukeysParam){
                    
                    tmpWgh*=pwr2(1.0-doubleVar/sqrdTukeysParam);
                    
                    MM[0]=tmpWgh*MM[0];
                    MM[1]=tmpWgh*MM[1];
                    MM[2]=tmpWgh*MM[2];
                    
                    modelm[0]+=MM[0];
                    modelm[1]+=MM[1];
                    modelm[2]+=MM[2];
                    
                    datam[0]+=tmpWgh*dataOriginalTemp[0];
                    datam[1]+=tmpWgh*dataOriginalTemp[1];
                    datam[2]+=tmpWgh*dataOriginalTemp[2];
                    
                    /* Compute the outer product */
                    Spx[0]+= MM[0]*dataOriginalTemp[0];
                    Spx[1]+= MM[1]*dataOriginalTemp[0];
                    Spx[2]+= MM[2]*dataOriginalTemp[0];
                    Spx[3]+= MM[0]*dataOriginalTemp[1];
                    Spx[4]+= MM[1]*dataOriginalTemp[1];
                    Spx[5]+= MM[2]*dataOriginalTemp[1];
                    Spx[6]+= MM[0]*dataOriginalTemp[2];
                    Spx[7]+= MM[1]*dataOriginalTemp[2];
                    Spx[8]+= MM[2]*dataOriginalTemp[2];
                    
                    quasum+=tmpWgh;
                    
                }
                
            }
            
        }
        
        if (quasum<1e-12){
            break;
        }
        
        /* Obtain the cross-covariance matrix, Spx */
        
        modelm[0]=modelm[0]/quasum;
        modelm[1]=modelm[1]/quasum;
        modelm[2]=modelm[2]/quasum;
        
        datam[0]=datam[0]/quasum;
        datam[1]=datam[1]/quasum;
        datam[2]=datam[2]/quasum;
        
        Spx[0]-=quasum*(modelm[0]*datam[0]);
        Spx[1]-=quasum*(modelm[1]*datam[0]);
        Spx[2]-=quasum*(modelm[2]*datam[0]);
        Spx[3]-=quasum*(modelm[0]*datam[1]);
        Spx[4]-=quasum*(modelm[1]*datam[1]);
        Spx[5]-=quasum*(modelm[2]*datam[1]);
        Spx[6]-=quasum*(modelm[0]*datam[2]);
        Spx[7]-=quasum*(modelm[1]*datam[2]);
        Spx[8]-=quasum*(modelm[2]*datam[2]);
        
        k=1;
        if (ii<1){
            
            doubleVar=0.0;
            for (j=0;j<9;j++){
                doubleVar+=Spx[j]*Spx[j];
            }
            
            doubleVar=doubleVar/pwr2(quasum);
            
            /* Check if Spx is in too bad condition in the first ICP-iteration  */
            if (doubleVar<1e-4){
                k=0;
                translationVector[0]=modelm[0]-datam[0];
                translationVector[1]=modelm[1]-datam[1];
                translationVector[2]=modelm[2]-datam[2];
            }
        }
        
        if (k) {
            
            /* calculate R using [U,Sigma,V]=svd(Spx) */
            
            SpxTSpx[0]=Spx[0]*Spx[0]+Spx[1]*Spx[1]+Spx[2]*Spx[2];
            SpxTSpx[1]=Spx[3]*Spx[3]+Spx[4]*Spx[4]+Spx[5]*Spx[5];
            SpxTSpx[2]=Spx[6]*Spx[6]+Spx[7]*Spx[7]+Spx[8]*Spx[8];
            SpxTSpx[3]=Spx[0]*Spx[3]+Spx[1]*Spx[4]+Spx[5]*Spx[2];
            SpxTSpx[4]=Spx[3]*Spx[6]+Spx[4]*Spx[7]+Spx[5]*Spx[8];
            SpxTSpx[5]=Spx[0]*Spx[6]+Spx[1]*Spx[7]+Spx[2]*Spx[8];
            
            
            /*   CharacteristicPolynomial  sigma^3-A*sigma^2-B*sigma+C   */
            
            A=SpxTSpx[2]+SpxTSpx[1]+SpxTSpx[0];
            B=SpxTSpx[5]*SpxTSpx[5]+SpxTSpx[4]*SpxTSpx[4]-SpxTSpx[2]*SpxTSpx[1]-SpxTSpx[2]*SpxTSpx[0]+SpxTSpx[3]*SpxTSpx[3]-SpxTSpx[1]*SpxTSpx[0];
            C=-2*SpxTSpx[5]*SpxTSpx[3]*SpxTSpx[4]+SpxTSpx[5]*SpxTSpx[5]*SpxTSpx[1]+SpxTSpx[4]*SpxTSpx[4]*SpxTSpx[0]+SpxTSpx[2]*SpxTSpx[3]*SpxTSpx[3]-SpxTSpx[2]*SpxTSpx[1]*SpxTSpx[0];
            
            sqrtA23B=sqrt(A*A+3*B);
            
            x0=(A-sqrtA23B)/3.0;
            f0=(x0*x0-A*x0-B)*x0;
            
            SIGMA[2]=MIN(MAX((x0*(C+2*f0-2*sqrt(f0*(f0+C)))/f0), 0), 0.5*x0);
            
            x0=(A+sqrtA23B)/3.0;
            f0=x0*x0*x0-A*x0*x0-B*x0+C;
            
            SIGMA[0]=MAX(MIN((x0*(B*A-C)+2*x0*f0+2*(x0-A)*sqrt(f0*(f0+B*A-C))-A*f0)/(f0+B*A-C), A), 0.5*(A+x0));
            
            for (k=0;k<3;k++){
                
                if (k==0){
                    j=0;
                }
                else if (k==1){
                    j=2;
                }
                else if (k==2){
                    j=1;
                    SIGMA[1]=A-SIGMA[0]-SIGMA[2];
                }
                
                /* Newton-Raphson */
                
                for (i=0;i<50;i++){
                    SIp=SIGMA[j];
                    difl=(-3*SIGMA[j]+2*A)*SIGMA[j]+B;
                    if (fabs(difl)>1e-15){
                        SIGMA[j]=((-2*SIGMA[j]+A)*SIGMA[j]*SIGMA[j]+C)/difl;
                        if (fabs(SIGMA[j]-SIp)<1e-25){
                            break;
                        }
                    }
                    else {
                        break;
                    }
                }
            }
            
            k=0;
            if (fabs(SIGMA[1]-SIGMA[0])<1e-12){
                k=1;
            }
            
            /* eigenvalues found, corresponding eigenvectors V[i] ... */
            
            for (i=0;i<3;i++){
                
                invmat[0]=SpxTSpx[2]*SpxTSpx[1]-SpxTSpx[1]*SIGMA[i]-SIGMA[i]*SpxTSpx[2]+SIGMA[i]*SIGMA[i]-SpxTSpx[4]*SpxTSpx[4];
                invmat[1]=SpxTSpx[4]*SpxTSpx[5]-SpxTSpx[3]*SpxTSpx[2]+SpxTSpx[3]*SIGMA[i];
                invmat[2]=SpxTSpx[3]*SpxTSpx[4]-SpxTSpx[5]*SpxTSpx[1]+SpxTSpx[5]*SIGMA[i];
                invmat[3]=SpxTSpx[2]*SpxTSpx[0]-SpxTSpx[0]*SIGMA[i]-SIGMA[i]*SpxTSpx[2]+SIGMA[i]*SIGMA[i]-SpxTSpx[5]*SpxTSpx[5];
                invmat[4]=-SpxTSpx[4]*SpxTSpx[0]+SpxTSpx[4]*SIGMA[i]+SpxTSpx[3]*SpxTSpx[5];
                invmat[5]=SpxTSpx[1]*SpxTSpx[0]-SpxTSpx[0]*SIGMA[i]-SpxTSpx[1]*SIGMA[i]+SIGMA[i]*SIGMA[i]-SpxTSpx[3]*SpxTSpx[3];
                
                if (i<2){
                    V[3*i]=invmat[0];
                    V[3*i+1]=invmat[1];
                    V[3*i+2]=invmat[2];
                    
                    if (k){
                        if (i==1){
                            V[3]=invmat[1];
                            V[4]=invmat[3];
                            V[5]=invmat[4];
                            
                            doubleVar=V[3]*V[0]+V[4]*V[1]+V[5]*V[2];
                            V[3]=V[3]-doubleVar*V[0];
                            V[4]=V[4]-doubleVar*V[1];
                            V[5]=V[5]-doubleVar*V[2];
                        }
                    }
                    
                }
                else {
                    
                    /* Eigen vectors corresponding to symmetric positiv definite matrices
                     * are orthogonal. */
                    
                    V[6]=V[1]*V[5]-V[2]*V[4];
                    V[7]=V[2]*V[3]-V[0]*V[5];
                    V[8]=V[0]*V[4]-V[1]*V[3];
                }
                
                for (j=0;j<10;j++){
                    
                    MM[0]=V[3*i];
                    MM[1]=V[3*i+1];
                    MM[2]=V[3*i+2];
                    
                    V[3*i]=invmat[0]*MM[0]+invmat[1]*MM[1]+invmat[2]*MM[2];
                    V[3*i+1]=invmat[1]*MM[0]+invmat[3]*MM[1]+invmat[4]*MM[2];
                    V[3*i+2]=invmat[2]*MM[0]+invmat[4]*MM[1]+invmat[5]*MM[2];
                    
                    if (k){
                        if (i==1){
                            doubleVar=V[3]*V[0]+V[4]*V[1]+V[5]*V[2];
                            V[3]=V[3]-doubleVar*V[0];
                            V[4]=V[4]-doubleVar*V[1];
                            V[5]=V[5]-doubleVar*V[2];
                        }
                    }
                    
                    doubleVar=sqrt(pwr2(V[3*i])+pwr2(V[3*i+1])+pwr2(V[3*i+2]));
                    
                    V[3*i]=V[3*i]/doubleVar;
                    V[3*i+1]=V[3*i+1]/doubleVar;
                    V[3*i+2]=V[3*i+2]/doubleVar;
                    
                    if (j>2){
                        if ((pwr2(V[3*i]-MM[0])+pwr2(V[3*i+1]-MM[1])+pwr2(V[3*i+2]-MM[2]))<1e-29){
                            break;
                        }
                    }
                }
            }
            
            /* singular values & V[i] of Spx found, U[i] ... */
            
            SpxTSpx[0]=Spx[0]*Spx[0]+Spx[3]*Spx[3]+Spx[6]*Spx[6];
            SpxTSpx[1]=Spx[1]*Spx[1]+Spx[4]*Spx[4]+Spx[7]*Spx[7];
            SpxTSpx[2]=Spx[2]*Spx[2]+Spx[5]*Spx[5]+Spx[8]*Spx[8];
            SpxTSpx[3]=Spx[0]*Spx[1]+Spx[3]*Spx[4]+Spx[6]*Spx[7];
            SpxTSpx[4]=Spx[1]*Spx[2]+Spx[4]*Spx[5]+Spx[7]*Spx[8];
            SpxTSpx[5]=Spx[0]*Spx[2]+Spx[3]*Spx[5]+Spx[6]*Spx[8];
            
            for (i=0;i<3;i++){
                
                invmat[0]=SpxTSpx[2]*SpxTSpx[1]-SpxTSpx[1]*SIGMA[i]-SIGMA[i]*SpxTSpx[2]+SIGMA[i]*SIGMA[i]-SpxTSpx[4]*SpxTSpx[4];
                invmat[1]=SpxTSpx[4]*SpxTSpx[5]-SpxTSpx[3]*SpxTSpx[2]+SpxTSpx[3]*SIGMA[i];
                invmat[2]=SpxTSpx[3]*SpxTSpx[4]-SpxTSpx[5]*SpxTSpx[1]+SpxTSpx[5]*SIGMA[i];
                invmat[3]=SpxTSpx[2]*SpxTSpx[0]-SpxTSpx[0]*SIGMA[i]-SIGMA[i]*SpxTSpx[2]+SIGMA[i]*SIGMA[i]-SpxTSpx[5]*SpxTSpx[5];
                invmat[4]=-SpxTSpx[4]*SpxTSpx[0]+SpxTSpx[4]*SIGMA[i]+SpxTSpx[3]*SpxTSpx[5];
                invmat[5]=SpxTSpx[1]*SpxTSpx[0]-SpxTSpx[0]*SIGMA[i]-SpxTSpx[1]*SIGMA[i]+SIGMA[i]*SIGMA[i]-SpxTSpx[3]*SpxTSpx[3];
                
                if (i<2){
                    U[3*i]=invmat[0];
                    U[3*i+1]=invmat[1];
                    U[3*i+2]=invmat[2];
                    
                    if (k){
                        if (i==1){
                            U[3]=invmat[1];
                            U[4]=invmat[3];
                            U[5]=invmat[4];
                            
                            doubleVar=U[3]*U[0]+U[4]*U[1]+U[5]*U[2];
                            U[3]=U[3]-doubleVar*U[0];
                            U[4]=U[4]-doubleVar*U[1];
                            U[5]=U[5]-doubleVar*U[2];
                        }
                    }
                    
                }
                else {
                    
                    /* Eigen vectors corresponding to symmetric positiv definite matrices
                     * are orthogonal. */
                    
                    U[6]=U[1]*U[5]-U[2]*U[4];
                    U[7]=U[2]*U[3]-U[0]*U[5];
                    U[8]=U[0]*U[4]-U[1]*U[3];
                }
                
                for (j=0;j<10;j++){
                    
                    MM[0]=U[3*i];
                    MM[1]=U[3*i+1];
                    MM[2]=U[3*i+2];
                    
                    U[3*i]=invmat[0]*MM[0]+invmat[1]*MM[1]+invmat[2]*MM[2];
                    U[3*i+1]=invmat[1]*MM[0]+invmat[3]*MM[1]+invmat[4]*MM[2];
                    U[3*i+2]=invmat[2]*MM[0]+invmat[4]*MM[1]+invmat[5]*MM[2];
                    
                    if (k){
                        if (i==1){
                            doubleVar=U[3]*U[0]+U[4]*U[1]+U[5]*U[2];
                            U[3]=U[3]-doubleVar*U[0];
                            U[4]=U[4]-doubleVar*U[1];
                            U[5]=U[5]-doubleVar*U[2];
                        }
                    }
                    
                    doubleVar=sqrt(pwr2(U[3*i])+pwr2(U[3*i+1])+pwr2(U[3*i+2]));
                    
                    U[3*i]=U[3*i]/doubleVar;
                    U[3*i+1]=U[3*i+1]/doubleVar;
                    U[3*i+2]=U[3*i+2]/doubleVar;
                    
                    if (j>2){
                        if ((pwr2(U[3*i]-MM[0])+pwr2(U[3*i+1]-MM[1])+pwr2(U[3*i+2]-MM[2]))<1e-29){
                            break;
                        }
                    }
                    
                }
                
            }
            
            k=0;
            for (i=0;i<9;i+=3){
                A=(Spx[0]*V[i]+Spx[3]*V[i+1]+Spx[6]*V[i+2])*U[i]+(Spx[1]*V[i]+Spx[4]*V[i+1]+Spx[7]*V[i+2])*U[i+1]+(Spx[2]*V[i]+Spx[5]*V[i+1]+Spx[8]*V[i+2])*U[i+2];
                if (A<0){
                    k++;
                    U[i]=-U[i];
                    U[i+1]=-U[i+1];
                    U[i+2]=-U[i+2];
                }
            }
            
            /* Get R=U*diag([1,1,det(U*V')])*V' */
            
            if (k==0 || k==2){           /* det(U*V')=+1 */
                for (i=0;i<3;i++){
                    for (j=0;j<3;j++){
                        rotationMatrix[i+3*j]=U[i]*V[j]+U[i+3]*V[j+3]+U[i+6]*V[j+6];
                    }
                }
            }
            else{                        /* det(U*V')=-1 */
                for (i=0;i<3;i++){
                    for (j=0;j<3;j++){
                        rotationMatrix[i+3*j]=U[i]*V[j]+U[i+3]*V[j+3]-U[i+6]*V[j+6];
                    }
                }
            }
            
            /* Get T=modelm-R*datam */
            
            translationVector[0]=modelm[0]-rotationMatrix[0]*datam[0]-rotationMatrix[3]*datam[1]-rotationMatrix[6]*datam[2];
            translationVector[1]=modelm[1]-rotationMatrix[1]*datam[0]-rotationMatrix[4]*datam[1]-rotationMatrix[7]*datam[2];
            translationVector[2]=modelm[2]-rotationMatrix[2]*datam[0]-rotationMatrix[5]*datam[1]-rotationMatrix[8]*datam[2];
            
        }
        
        if (bol){
            numRndInIters+=rndIncrease;
            if (numRndInIters>=numDataPoints){
                bol=0;
            }
        }
        
    }
    
    free(Nvertex_1);
    free(NvertexSq);
    free(Nvertex_1xn);
    free(Nvertex_1xnxn);
    free(dNvertex);
    free(dNvertex_1);
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {
    
    if (nrhs != 13){
        mexErrMsgTxt("13 input arguments are required");
    } else if (nlhs != 2) {
        mexErrMsgTxt("2 output arguments are required");
    }
    
    rotationMatrixOut = mxCreateDoubleMatrix(3, 3, mxREAL);
    transVecOut = mxCreateDoubleMatrix(3, 1, mxREAL);
    
    icpDVG(mxGetPr(data), mxGetPr(quality), (unsigned int)mxGetN(data), (unsigned int*)mxGetPr(randvec), (unsigned int)mxGetScalar(sizerand), (unsigned int)mxGetScalar(randincrease), (unsigned int)mxGetScalar(maxiter), mxGetScalar(kTukeySq), mxGetPr(model), mxGetPr(dirvec), mxGetPr(normal), mxGetPr(pdir), mxGetPr(sdir), DVGdata, mxGetPr(rotationMatrixOut), mxGetPr(transVecOut));
    
    return;
    
}
