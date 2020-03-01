#include "ukf.h"
#include <math.h>
#include <stdlib.h>
#if KAL_LAYOUT == lahColMajor 
#include <string.h>
#endif

/*To Do: Change Kalman struct to distinguish between simple
 * and unscented
 * Unscented needs: Paramters alpha beta kappa
 * Simple needs: Transition function with alpha*/
lah_Return ukal_correct(ukal_Filter *ukal, const lah_mat *z)
{
    lah_Return result = lahReturnOk;
    
    lah_index j;
    lah_value temp;
    
    lah_mat *Measurements = ukal->Measurements;
    lah_mat *measurement = ukal->measurement;
    lah_mat *sigmaPoint = ukal->sigmaPoint;
    lah_mat *SigmaPoints = ukal->SigmaPoints;
    
    lah_mat *CrossCovariance = ukal->CrossCovariance;
    
    /*--- Apply measurement model to all sigma points ---*/
    /* At the same time, we calculate the mean z = sum(Weights * meas)*/
    result += ukal->h(ukal->meas_zero, ukal->x, ukal->optData);
    lah_matAdd(ukal->z, 0.0, ukal->z, ukal->weight_mean, ukal->meas_zero);
    
    sigmaPoint->data = SigmaPoints->data; /* Set on first col*/
    measurement->data = Measurements->data; /* Reset col slice */
    for (j = 0; j < Measurements->nC; j++)
    {
        result += ukal->h(measurement, sigmaPoint, ukal->optData);
        /* Update predicted measurement */
        lah_matAdd(ukal->z, 1.0 , ukal->z, ukal->weight, measurement);
        
        measurement->data += Measurements->incCol;
        sigmaPoint->data += SigmaPoints->incCol;
    }
    if (result != lahReturnOk)
        return lahReturnFPError;
    
    /*--- Substract mean of z from all Measurements ---*/
    lah_matAdd(ukal->meas_zero, 1.0, ukal->meas_zero, -1.0, ukal->z);
    measurement->data = Measurements->data; /* Reset col slice */
    for (j = 0; j < Measurements->nC ; j++)
    {
        /* lah_matAdd(measurement, 1.0, measurement, -1.0, ukal->z); */
        lah_matAdd(measurement, sqrt(ukal->weight), measurement, -1.0 * sqrt(ukal->weight), ukal->z);
        measurement->data += Measurements->incCol;
    }
    
    /*--- Create the CrossCovariance matrix ---*/
    /* Cov(x, z) = sum_j (W_j (x_j -x)(z_j -z)^t) */
    
    /* Treat first iteration seperately to initiliaze CrossCovariance*/
    lah_matMul(lahNorm, lahTrans, 0.0, ukal->weight_cov, ukal->CrossCovariance, ukal->sigma_zero, ukal->meas_zero);
    
    /* Reset col slices */
    measurement->data = Measurements->data; 
    sigmaPoint->data = SigmaPoints->data;
    for (j = 0; j < Measurements->nC; j++)
    {
        /*lah_matMul(lahNorm, lahTrans, 1.0, ukal->weight, ukal->CrossCovariance, sigmaPoint, measurement);*/
        lah_matMul(lahNorm, lahTrans, 1.0, 0.0, ukal->CrossCovariance, sigmaPoint, measurement);
        measurement->data += Measurements->incCol;
        sigmaPoint->data += SigmaPoints->incCol;
    }
    
    /*--- QR factorization of Augmented Measurment Matrix ---*/
    /* The left part of the matrix holds the 2*nState 
     * virtual measurements (without meas_zero). 
     * The right part holds Rsqrt  */
    
    /*Copy Rsqrt to right part of Matrix, 
     * will be overwritten by RootCovariance */
    
    #if LAH_LAYOUT == lahColMajor 
    memcpy(ukal->RootCovariance_meas->data, ukal->Rsqrt->data, 
           measurement->nR * measurement->nR * sizeof(lah_value));
    #else   
    for (i = 0; i < ukal->RootCovariance_meas->nR; i++)
    {
        for(j = 0; j < ukal->RootCovariance_meas->nC; j++)
        {
            LAH_ENTRY(ukal->RootCovariance_meas, i, j) 
                = LAH_ENTRY(ukal->Rsqrt, i, j);
        }
    }
    #endif
    
    lah_QR(ukal->AugmentedMeas_trans, ukal->QR_tau);
    
    /*--- Do a rank-one Cholesky update for zeroth virt. measurement ---*/
    if (ukal->weight_cov < 0)
        temp = -1.0 * sqrt(fabs(ukal->weight_cov));
    else
        temp = sqrt(ukal->weight_cov);
    
    lah_cholUpdate(ukal->RootCovariance_meas, ukal->meas_zero, temp);
    
    /*--------------------*/
    /* Correction step    */
    /*--------------------*/
    
    /*--- Calculate Kalman gain ---*/    
    /* (RootCov'\(RootCov\CrossCov'))' 
     * equivalent to (CrossCov/RootCov')/ RootCov 
     */
    lah_forwardSub(ukal->CrossCov_trans, ukal->RootCovariance_meas);
    /*After this operation CrossCovariance will hold the Gain Matrix*/
    lah_forwardSub(ukal->CrossCov_trans, ukal->RootCov_trans);

    /*--- Apply Update ---*/
    /* Update mean: x = x + Gain * (z - zp) */
    /* Use meas_zero as workspace  for error: meas_zero = z - zp */
    lah_matAdd(ukal->meas_zero, 1.0, z, -1.0, ukal->z);
    lah_matMul(lahNorm, lahNorm, 1.0, 1.0, ukal->x, CrossCovariance, ukal->meas_zero);

    /*--- Calculate Cholesky update matrix ---*/
    /* Use Gain as workspace and save result there*/
    lah_matMul(lahNorm, lahNorm, 0.0, 1.0, ukal->Gain, 
               CrossCovariance, ukal->RootCovariance_meas);

    /*--- Update Rootcovariance using nMeas rank-one Cholesky downdates ---*/
    sigmaPoint->data = ukal->Gain->data; 
    for(j = 0; j < ukal->Gain->nC; j++) 
    {
        result += lah_cholUpdate(ukal->RootCovariance_meas, sigmaPoint, -1.0);
        sigmaPoint->data += ukal->Gain->incCol;
    }
    return result;
}
