#ifndef _UKF_H_
#define _UKF_H_
#include <lah.h>

#ifdef __cplusplus
extern "C" {
#endif
/*****************************************************************************/
/*  Implementation of the unscented Kalman Filter                            */
/*****************************************************************************/

/*--- Definitions for unscented Kalman filter -------------------------------*/

/* ukal_StateFun
 * -------------
 * Measurement function for the unscented Kalman filter.
 * Describes the function z = h(x), relating state to measurement.
 */
typedef lah_Return (*ukal_StateFun)(lah_mat *xp, const lah_mat *x,
                                   const lah_mat *u, void *optData);

/* ukal_MeasFun
 * -------------
 * Measurement function for the unscented Kalman filter.
 * Describes the function z = h(x), relating state to measurement.
 */
typedef lah_Return (*ukal_MeasFun)(lah_mat *zp, const lah_mat *x, 
                                   void *optData);

/* ukal_Filter: Data struct to compute the unscented Kalman filter          */
/* Uses the square root implementation, where not the Covariance, but the   */
/* SquareRoot of the Covariance is updated                                  */
typedef struct
{
    lah_mat *x;              /* state vector mean: nState x 1 */
    lah_mat *xwork;          /* workspace for mean, also used for sigma_0: 
                                nState x 1 */

    lah_mat *z;              /* measurement vector mean: nMeas x 1 */
    
    /* Propagation */
    ukal_StateFun f;         /* state prediction function */
	ukal_MeasFun h;	         /* measurement prediction function */
    void *optData;           /* optional user data for propagation */
                               
    /* Matrices */
    lah_mat *AugmentedSigma; /* Matrix holding data for sigmaPoints and 
                                RootCovariance_state*/
    lah_mat *AugmentedMeas;  /* Matrix holding data for virtual measurements
                                and RootCovariance_meas*/
    
    lah_mat *Qsqrt;          /* Sqrt of process noise covariance */
    lah_mat *Rsqrt;          /* Sqrt of measurement noise covariance */  
    
    lah_mat *sigma_zero;     /* zeroth sigmaPoint, treated seperately: nState x 1 */      
    lah_mat *meas_zero;      /* zeroth virtual meas, treated seperately: nMeas x 1 */      
    
    lah_mat *Gain;           /* nState x nMeas */
    lah_mat *CrossCovariance; /* nState * nMeas */
    
    /* Array for QR */
    lah_value *QR_tau;
    
    /* Matrices without own data */
    lah_mat *RootCovariance_state; /* nState * nState */
    lah_mat *RootCovariance_meas; /* nMeas * nMeas */
    
    lah_mat *HelperMatrix;     /* nState x nState Matrix used 
                                  in redrawSigmaPoints */
    
    /* Sigma Points */
    lah_mat *SigmaPoints;    /* nState x 2*nState */
    lah_mat *sigmaPoint;     /* nState x 1 */
    lah_mat *Measurements;   /* nMeas x 2*nState */
    lah_mat *measurement;    /* nMeas x 1 */
    
    /*Transposed Views*/
    lah_mat *CrossCov_trans;
    lah_mat *RootCov_trans;
    lah_mat *Gain_trans;
    lah_mat *AugmentedMeas_trans;
    lah_mat *AugmentedSigma_trans; 
    
    /* Scaling Weights */
    lah_value weight;      /* Scaling Weights for SigmaPoints */
    lah_value weight_mean;  /* Weight for means */
    lah_value weight_cov;   /* Weight for covariance */
} ukal_Filter;

/*--------------------------------------------------------------------------*/

/*--- Unscented Kalman-Filter functions ------------------------------------*/ 
/* Square-root implementation of the unscented Kalman-Filter                */
/* For details see "Rudolph van der Merwe and Eric A. Wan.                  */ 
/* - The square-root unscented Kalman filter for state and                  */
/* parameter-estimation"                                                    */
/* This implementation uses scaling weights for sigma Points, which are     */
/* defined by three parameters alpha, beta and kappa in ukal_create         */
/*--------------------------------------------------------------------------*/

/* ukal_predict:
 * ------------
 * This function predicts the next filter state using the current state, 
 * the given input u and the given process covariance. It will call 
 * the function f internally to infer the process state.
 */
lah_Return ukal_predict(ukal_Filter *ukal, const lah_mat *u);

/* ukal_predict:
 * ------------
 * This function corrects the current filter state by using the given real measurement z and
 * the given measurement process noise covariance R.
 */
lah_Return ukal_correct(ukal_Filter *ukal, const lah_mat *z);


lah_Return ukal_redrawSigmaPoints(ukal_Filter *ukal);

/*--- Utility/Helper functions ---------------------------------------------*/

/* Allocates all data for the unscented */
ukal_Filter *ukal_create(lah_index nState, lah_index nMeas,
                         ukal_StateFun f, ukal_MeasFun h, 
                         void *optData);

/* Allocates all data for the unscented */
lah_Return ukal_init(ukal_Filter *ukal, lah_mat *Q, lah_mat *R,
                     lah_mat *x_0, lah_value mean_value,
                     lah_value alpha, lah_value beta, lah_value kappa);

/* Allocates all data for the unscented */
ukal_Filter *ukal_free(ukal_Filter *ukal);

/*--------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* _UKF_H_ */
