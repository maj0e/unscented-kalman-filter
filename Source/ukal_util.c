#include "ukf.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include <stdio.h>

/* allocation for the feature vector feat */
ukal_Filter *ukal_create(lah_index nState, lah_index nMeas,
                         ukal_StateFun f, ukal_MeasFun h, void *optData)
{  
    ukal_Filter *ukal;
    
    /* Check parameters*/
    if (f == NULL || h == NULL) 
        return NULL;
        
    /* Allocate the struct */
    ukal = calloc (1, sizeof (ukal_Filter));
    if (!ukal) return (NULL);     /* out of memory */
    
    /*--- Allocate matrices of the Kalman Filter ---*/
    /* State vector*/
    ukal->x = lah_matAlloc(1, nState, 1);
    ukal->xwork = lah_matAlloc(1, nState, 1);
    
    /* Pred. meas. vector*/
    ukal->z  = lah_matAlloc(1, nMeas, 1);
    
    /*--- Set Propagation functions ---*/
    ukal->f = f;
    ukal->h = h;
    ukal->optData = optData;
    
    /*--- Allocate Matrices ---*/
    /* Matrix holding data for sigmaPoints and RootCovariance_state*/
    ukal->AugmentedSigma = lah_matAlloc(3 * nState, nState, 1); 
    /* Matrix holding data for virtual measurements and RootCovariance_meas*/
    ukal->AugmentedMeas = lah_matAlloc(2 * nState + nMeas, nMeas, 1);  
    
    ukal->Qsqrt = lah_matAlloc(nState, nState, 1); /* Sqrt of process noise covariance */
    ukal->Rsqrt = lah_matAlloc(nMeas, nMeas, 1); /* Sqrt of measurement noise covariance */  
    
    ukal->sigma_zero = lah_matAlloc(1, nState, 1);     /* zeroth sigmaPoint, treated seperately: nState x 1 */      
    ukal->meas_zero = lah_matAlloc(1, nMeas, 1);      /* zeroth virtual meas, treated seperately: nMeas x 1 */      
    
    ukal->Gain = lah_matAlloc(nMeas, nState, 1);           /* nState x nMeas */
    ukal->CrossCovariance = lah_matAlloc(nMeas, nState, 1); /* nState * nMeas */
    
    /* QR array */
    if (nMeas > nState)
        ukal->QR_tau = malloc(nMeas * sizeof(lah_value));
    else
        ukal->QR_tau = malloc(nState * sizeof(lah_value));
    
    /*--- Alllocate Matrix Views ---*/
    ukal->RootCovariance_state = lah_matAlloc(nState, nState, 0); /* nState * nState */
    ukal->RootCovariance_meas = lah_matAlloc(nMeas, nMeas, 0);   /* nMeas * nMeas */
    
    ukal->HelperMatrix = lah_matAlloc(nState, nState, 0);
    
    /* Sigma Points */
    ukal->SigmaPoints = lah_matAlloc(2 * nState, nState, 0);    /* nState x 2*nState */
    ukal->sigmaPoint = lah_matAlloc(1, nState, 0);     /* nState x 1 */
    ukal->Measurements = lah_matAlloc(2 * nState, nMeas, 0);   /* nMeas x 2*nState */
    ukal->measurement = lah_matAlloc(1, nMeas, 0);    /* nMeas x 1 */
    
    /*--- Configure Matrix Views ---*/
    ukal->RootCovariance_state->data 
        = ukal->AugmentedSigma->data + 2 * nState * ukal->AugmentedSigma->incCol;
    ukal->RootCovariance_state->incRow = ukal->AugmentedSigma->incRow;
    ukal->RootCovariance_state->incCol = ukal->AugmentedSigma->incCol;
    
    ukal->RootCovariance_meas->data
        = ukal->AugmentedMeas->data + 2 * nState * ukal->AugmentedMeas->incCol;
    ukal->RootCovariance_meas->incRow = ukal->AugmentedMeas->incRow;
    ukal->RootCovariance_meas->incCol = ukal->AugmentedMeas->incCol;
    
    ukal->SigmaPoints->data = ukal->AugmentedSigma->data;
    ukal->SigmaPoints->incRow = ukal->AugmentedSigma->incRow;
    ukal->SigmaPoints->incCol = ukal->AugmentedSigma->incCol;
    
    ukal->Measurements->data = ukal->AugmentedMeas->data;
    ukal->Measurements->incRow = ukal->AugmentedMeas->incRow;
    ukal->Measurements->incCol = ukal->AugmentedMeas->incCol;
    
    /* Create tranposed Views */
    ukal->CrossCov_trans = lah_matTrans(ukal->CrossCovariance);
    ukal->RootCov_trans = lah_matTrans(ukal->RootCovariance_meas);
    ukal->Gain_trans = lah_matTrans(ukal->Gain);
    ukal->AugmentedMeas_trans = lah_matTrans(ukal->AugmentedMeas);
    ukal->AugmentedSigma_trans = lah_matTrans(ukal->AugmentedSigma); 
    
    /* Configure column slices */
    ukal->measurement->incRow = ukal->Measurements->incRow;
    ukal->measurement->incCol = 0;                      /* Needed for repeated matrix view */
    ukal->measurement->data = ukal->Measurements->data;
    ukal->sigmaPoint->incRow = ukal->SigmaPoints->incRow;
    ukal->sigmaPoint->incCol = 0;                       /* Needed for repeated matrix view */
    ukal->sigmaPoint->data = ukal->SigmaPoints->data;

    /* Do checks if everything worked*/
    if (!ukal->x || !ukal->xwork || !ukal->z 
        || !ukal->AugmentedSigma || !ukal->AugmentedMeas
        || !ukal->Qsqrt || !ukal->Rsqrt
        || !ukal->sigma_zero || !ukal->meas_zero
        || !ukal->Gain || !ukal->CrossCovariance
        || !ukal->RootCovariance_state || !ukal->RootCovariance_meas
        || !ukal->HelperMatrix || !ukal->CrossCov_trans
        || !ukal->RootCov_trans || !ukal->Gain_trans    
        || !ukal->SigmaPoints || !ukal->sigmaPoint 
        || !ukal->Measurements || !ukal->measurement )
    {
        return ukal_free(ukal);
    }
    else
        return ukal;
}

ukal_Filter *ukal_free(ukal_Filter *ukal)
{
    if (!ukal) return NULL;      /* do nothing if kal already NULL */
    
    /* Free all data of the unscented Kalman-Filter 
     * According to C standard free(ptr) does nothing 
     * if ptr is already NULL
     */
    /* Free Matrices */
    lah_matFree(ukal->x);
    lah_matFree(ukal->xwork);
    lah_matFree(ukal->z);
    lah_matFree(ukal->AugmentedSigma);
    lah_matFree(ukal->AugmentedMeas);
    lah_matFree(ukal->Qsqrt);
    lah_matFree(ukal->Rsqrt);
    lah_matFree(ukal->sigma_zero);
    lah_matFree(ukal->meas_zero);
    lah_matFree(ukal->Gain);
    lah_matFree(ukal->CrossCovariance);
    
    free(ukal->QR_tau); 
    
    /* Free Matrix Views */
    free(ukal->RootCovariance_state);
    free(ukal->RootCovariance_meas);
    free(ukal->HelperMatrix);
    free(ukal->SigmaPoints);
    free(ukal->sigmaPoint);
    free(ukal->Measurements);
    free(ukal->measurement);
    free(ukal->CrossCov_trans);
    free(ukal->RootCov_trans);
    free(ukal->Gain_trans); 
    free(ukal->AugmentedMeas_trans);
    free(ukal->AugmentedSigma_trans); 
    
    /* free ukal itself*/
    free(ukal);
    return NULL; 
}

lah_Return ukal_init(ukal_Filter *ukal, lah_mat * Q, lah_mat *R,
                     lah_mat *x_0, lah_value mean_value,
                     lah_value alpha, lah_value beta, lah_value kappa)
{
    lah_Return result;
    lah_index i;
    lah_value scaling;
    lah_index nState = ukal->x->nR;    
    
    if (Q == NULL || R == NULL) 
        return lahReturnParameterError;
    if (Q->nR != nState || Q->nC != nState 
        || R->nR != ukal->z->nR || R->nC != ukal->z->nR) 
        return lahReturnParameterError;
    
    /* Compute the sqrt of R and Q: 
     * Use RootCovariance as workspace */
    printf("Init Sqrt error matrices\n");
    lah_matSqrt(Q, ukal->Qsqrt, ukal->RootCovariance_state);
    lah_matSqrt(R, ukal->Rsqrt, ukal->RootCovariance_meas);
    
    /*--- Compute weights ---*/
    printf("Compute Weights\n");
    scaling = alpha * alpha * (nState + kappa);
    printf("Scaling = %g\n", scaling);
    ukal->weight_mean = (scaling - nState)/scaling;
    printf("weight_mean = %g\n", ukal->weight_mean);
    ukal->weight_cov = ukal->weight_mean + (1.0 - (alpha*alpha) + beta);
    printf("weight_cov = %g\n", ukal->weight_cov);
    ukal->weight = 0.5 / scaling;
    printf("weight = %g\n", ukal->weight);
    
    /* Set x to initial x_0*/
    printf("Copy Initial x\n");
    memcpy(ukal->x->data, x_0, nState * sizeof(lah_value));
    
    /*Use x_0 now as workspace for x - mean_value*/
    for (i = 0; i < x_0->nR; i++)
        LAH_ENTRY(x_0, i, 0) = LAH_ENTRY(x_0, i, 0) - mean_value;
    
    /* Initialize Square Root matrix */
    /* RootCovariance = chol(E[(x - mean_value)(x - mean_value)^T])*/
    lah_matMul(lahNorm, lahTrans, 0.0, 1.0, ukal->RootCovariance_state, x_0, x_0);
    printf("Cholesky\n");
    result = lah_chol(ukal->RootCovariance_state, 1);
    /* Calculate sigma points and weigths */
    ukal_redrawSigmaPoints(ukal);
    
    return result;
}
