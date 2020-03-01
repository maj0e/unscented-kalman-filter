/*******************************************************************/
/**  Example to test the Kalman-Filter functions                  **/
/*******************************************************************/
#include "kalman.h"
#include "lah.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*the state prediction function: linear case for simplicity */

const lah_value accel   = 0.3;
const lah_value deltaT  = 0.1;
const lah_value sigma_x = 1;
const lah_value sigma_z = 3;

/* TODO: Set reasonable values for parameters */
#define ALPHA 0.2
#define BETA 2.0
#define KAPPA 0

struct optData 
{
    lah_index inits;
    lah_index initm;
    lah_mat *A;
    lah_mat* B; /* Control matrix*/
};

lah_Return timeEvolve(lah_mat* xp, const lah_mat *x,
                      const lah_mat *u, struct optData *Data)
{
	/* Set the state transition Matrix A
     * if not done yet*/
    lah_mat *B = Data->B;
    lah_mat *A = Data->A;
    lah_Return result;
    if (Data->inits == 0 )
    {
        A->data[0] = 1;
        A->data[A->incRow] = 0;
        A->data[A->incCol] = deltaT;
        A->data[A->incCol + A->incRow] = 1;
        
        B->data[0] = deltaT * deltaT /2;
        B->data[1] = deltaT;
        Data->inits = 1;
    }

    /* compute predicted state */
    result = lah_matMul(lahNorm, lahNorm, 0.0, 1.0 , xp, A, x);
    result += lah_matMul(lahNorm, lahNorm, 1.0, 1.0, xp, B, u);
        
    if (result != 0)
        return lahReturnMathError;
    else
        return lahReturnOk;
}

/*  measurement prediction function */
lah_Return measure(lah_mat *zp, const lah_mat *x,
                   struct optData *Data)
{
    /* compute the measurement from state x */
    zp->data[0] = x->data[0];

    return lahReturnOk;
}

int main ()
{
    const lah_index nState = 2;
    const lah_index nMeas = 1;
    struct optData Data;
    
    lah_index i;
    lah_value vel = 0, pos = 0;
    lah_value noise;
    lah_mat *x_0 = lah_matAlloc(1, 2, 1);
    lah_mat * u = lah_matAlloc(1, 1, 1);
    lah_mat * z = lah_matAlloc(1, 1, 1);
    lah_mat * R = lah_matAlloc(1, 1, 1);
    lah_mat * Q = lah_matAlloc(2, 2, 1);
    ukal_Filter* ukal;
    LAH_ENTRY(Q, 0, 0) = pow(sigma_x, 2) * pow(deltaT, 4) / 4;
    /*LAH_ENTRY(Q, 1, 0) = pow(sigma_x, 2) * pow(deltaT, 3) / 2;
    LAH_ENTRY(Q, 0, 1) = pow(sigma_x, 2) * pow(deltaT, 3) / 2; */ 
    LAH_ENTRY(Q, 1, 1) = pow(sigma_x, 2) * pow(deltaT, 2);

    /*Initialize the struct for prediction */
    Data.initm = 0;
    Data.inits = 0;
    Data.B = lah_matAlloc(1, 2, 1);
    Data.A = lah_matAlloc(2, 2, 1);
    printf("Create Unscented Kalman Filter\n");
    ukal = ukal_create(nState, nMeas,
                             (ukal_StateFun) timeEvolve, 
                             (ukal_MeasFun) measure, (void*) &Data);
    
    /* Initialization */ 
    /* input and process noise variables */
    u->data[0] = accel;
    R->data[0] = sigma_z * sigma_z;
    LAH_SETTYPE(Q, lahTypeDiagonal);
    LAH_SETTYPE(R, lahTypeDiagonal);
    
    
    LAH_ENTRY(x_0, 0, 0) = 0.75; /* Initial x_0 */
    LAH_ENTRY(x_0, 1, 0) = 0.5;
    
    
    
    printf("Initialize Unscented Kalman Filter\n");
    ukal_init(ukal, Q, R, x_0, 0.0, ALPHA, BETA, KAPPA);
    LAH_ENTRY(ukal->RootCovariance_state, 1, 1) = 1;
    lah_matPrint(ukal->RootCovariance_state, 1);
    lah_matPrint(ukal->RootCovariance_meas, 1);
    LAH_ENTRY(ukal->Qsqrt, 0, 0) = sqrt( pow(sigma_x, 2) * pow(deltaT, 4) / 4);
    /*LAH_ENTRY(Q, 1, 0) = pow(sigma_x, 2) * pow(deltaT, 3) / 2;
    LAH_ENTRY(Q, 0, 1) = pow(sigma_x, 2) * pow(deltaT, 3) / 2; */ 
    LAH_ENTRY(ukal->Qsqrt, 1, 1) =sqrt ( pow(sigma_x, 2) * pow(deltaT, 2));
    lah_matPrint(ukal->Qsqrt, 1);
    lah_matPrint(ukal->Rsqrt, 1);
    lah_matFree(Q);
    lah_matFree(R);
    lah_matFree(x_0);
    /* initialize random number generator */
    srand(0);

    /* print header */
    printf("Start Loop...\n");
    printf("i pos vel realpos realvel z zp\n");
    ukal_predict(ukal, u); /* Without these lines --> NaNs , search for origin*/
    /* loop over time and demonstrate Kalman filter */
    for (i = 1; i < 500; i++)
    {
        /* virtual measurement */
        noise = lah_gaussNoise();
        pos = u->data[0] / 2.0 * pow(i * deltaT, 2.0);
        vel = u->data[0] * i * deltaT;
		z->data[0] = pos + noise * sigma_z;
        
        /* correct filter state */
        ukal_correct(ukal, z);
        /* print out */
        printf("%d %g %g %g %g %g %g\n", i, 
                ukal->x->data[0], ukal->x->data[1], pos, vel,
                z->data[0], ukal->z->data[0]);

		u->data[0] = accel; /* unneccesary ?? */				
        /* predict the next filter state */
    	ukal_predict(ukal, u);
    }
    lah_matFree(u);
    lah_matFree(z);
    lah_matFree(Data.B);
    lah_matFree(Data.A);
    ukal_free(ukal);

    return 0;
}
