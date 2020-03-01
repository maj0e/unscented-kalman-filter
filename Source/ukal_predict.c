#include "ukf.h"
#include <math.h>
#include <stdlib.h>
#if KAL_LAYOUT == lahColMajor 
#include <string.h>
#endif

lah_Return ukal_predict(ukal_Filter *ukal, const lah_mat *u)
{
    lah_Return result = lahReturnOk;
    lah_index j;
    lah_value temp;
    
    lah_mat *sigmaPoint = ukal->sigmaPoint;
    lah_mat *SigmaPoints = ukal->SigmaPoints;
    lah_mat *xwork = ukal->xwork;
    lah_index nState= ukal->x->nR;
    
    if (ukal == NULL)
    {
        return lahReturnParameterError;
    }
    if (ukal->f == NULL)
        return lahReturnFPError;
    
    /*--- Propagate all sigma points ---*/
    /* Treat zeroth sigmaPoint seperately */
    result += ukal->f(ukal->xwork, ukal->sigma_zero, u, ukal->optData);
    lah_matAdd(ukal->x, 0.0 , ukal->x, ukal->weight_mean, xwork);
    memcpy(ukal->sigma_zero->data, xwork->data, nState * sizeof(lah_value));
    
    sigmaPoint->data = SigmaPoints->data; /* Set on first col*/
    for (j = 0; j < SigmaPoints->nC; j++)
    {   
        /* Propagate sigma point */
        result += ukal->f(ukal->xwork, sigmaPoint, u, ukal->optData);
        
        /* Update predicted mean*/
        lah_matAdd(ukal->x, 1.0 , ukal->x, ukal->weight, xwork);
        
        memcpy(sigmaPoint->data, xwork->data, nState * sizeof(lah_value));
        sigmaPoint->data += SigmaPoints->incCol;
    }
    if (result != lahReturnOk)
        return lahReturnFPError;
    
    /*------ Substract mean from all SigmaPoints ---------------------------*/
    lah_matAdd(ukal->sigma_zero, 1.0, ukal->sigma_zero, -1.0, ukal->x);
    sigmaPoint->data = SigmaPoints->data; /* Reset col slice */
    for (j = 0; j < SigmaPoints->nC ; j++)
    {
        /*lah_matAdd(sigmaPoint, 1.0, sigmaPoint, -1.0, ukal->x);*/
        lah_matAdd(sigmaPoint, sqrt(ukal->weight), sigmaPoint, -1.0 * sqrt(ukal->weight), ukal->x);
        sigmaPoint->data += SigmaPoints->incCol;
    }
    /*--- QR factorization of Augmented Sigma Matrix ---*/
    /* The left part of the matrix holds the 2*nState SigmaPoints 
     * (without sigma_zero). The right part holds Qsqrt  */
    
    /*Copy Qsqrt to right part of Matrix, 
     * will be overwritten by RootCovariance */
    #if KAL_LAYOUT == lahColMajor 
        memcpy(ukal->RootCovariance_state->data, ukal->Qsqrt->data, 
               nState * nState * sizeof(lah_value));
    #else
    for (i = 0; i < ukal->RootCovariance_state->nR; i++)
    {
        for(j = 0; j < ukal->RootCovariance_state->nC; j++)
        {
            LAH_ENTRY(ukal->RootCovariance_state, i, j) 
                = LAH_ENTRY(ukal->Qsqrt, i, j);
        }
    }
    #endif
    
    result += lah_QR(ukal->AugmentedSigma_trans, ukal->QR_tau);
    if (result != lahReturnOk)
        return result;
    
    /*------ Rank-one Cholesky update of root covariance matrix ---*/
    /* using zeroth sigma point */
    if (ukal->weight_mean < 0)
        temp = -1.0 * sqrt(fabs(ukal->weight_cov));
    else
        temp = sqrt(ukal->weight_cov);
    
    result += lah_cholUpdate(ukal->RootCovariance_state, ukal->sigma_zero, temp);
    if (result != lahReturnOk)
        return result;

    /*--- Redraw sigma points using a priori root covariance ---*/
    return ukal_redrawSigmaPoints(ukal);
}
