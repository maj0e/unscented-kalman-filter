#include "ukf.h"
#include <math.h>
#include <string.h>


#if 0
lah_Return ukal_calcSigmaPoints(kal_unscentedFilter *ukal)
{
    lah_index i,j;
    lah_index nDim = x->nR;
    lah_value scaling = alpha*alpha *(nDim + kappa); /* Put in constructor of struct, can be saved for all calculations*/

    /* Need to change incRow incCol to get Submatrices as workspace*/
    /*if something breaks, look here first :) */
    lah_mat *Psqrt = lah_matAlloc(nDim, nDim, KAL_LAYOUT, 0);
    lah_mat *work = lah_matAlloc(nDim, nDim, KAL_LAYOUT, 0);
    lah_mat *SigmaPoints_col = lah_matAlloc(1, nDim, KAL_LAYOUT, 0);*/
    
    Psqrt->data = SigmaPoints->data + SigmaPoints->incCol; /*start with second column */
    Psqrt->incRow = SigmaPoints->incRow;
    Psqrt->incCol = SigmaPoints->incCol;

    work->data = SigmaPoints->data + SigmaPoints->incCol * (nDim+1); /* use the last nDim columns */
    work->incRow = SigmaPoints->incRow;
    work->incCol = SigmaPoints->incCol;

    /*Calculate sqrt(scaling * P) and save in SigmaPoints matrix */
    lah_matSqrt(P, Psqrt, work);

    /* First sigma point --> Mean x */
    for (j = 0; j < ukal->x->nR; j++)
        mean->data[mean->incRow * j] = x->data[x->incRow * j];

    /* SigmaPoints_meancol is only one column holding the mean vector x.
     * Change interface to look like a Matrix 
     * --> mean vector repeated nDim times  */
    SigmaPoints_meancol->incCol = 0;
    SigmaPoints_meancol->nC = nDim;
    
    /*Set Sigma Points nDim+1 to 2nDim+1 */
    /* SigmaPoints = x - sqrt(scaling) P_sqrt*/
    lah_matAdd(work, 1, SigmaPoints_meancol, (-1) * sqrt(scaling), P_sqrt);

    /* Set Sigma Points 1 to nDim*/
    /* SigmaPoints = x - sqrt(scaling) P_sqrt*/
    lah_matAdd(Psqrt, 1, SigmaPoints_meancol, sqrt(scaling), P_sqrt);

    /* Calculate Weights */
    weights->data[0] = (scaling - nDim)/(nDim + scaling);
    for (j = 1; j <= 2 * nDim ; j++)
        weights->data[j] = 0.5 / scaling;

    /*Clean Up*/
    free(SigmaPoints_meancol);
    free(work);
    free(Psqrt);
}
#endif

lah_Return ukal_redrawSigmaPoints(ukal_Filter *ukal)
{
    lah_Return result;
    lah_mat *mean = ukal->x;
    lah_index nState = mean->nR;
    lah_mat *sigmaPoint = ukal->sigmaPoint;
    lah_mat *SigmaPoints = ukal->SigmaPoints;
    lah_mat *HelperMatrix = ukal->HelperMatrix;

    /* First sigma point --> Mean x */
    /* TODO: Maybe we can jsut swap pointers here: is mean needed afterwards???*/
    memcpy(ukal->sigma_zero->data, mean->data, nState * sizeof(lah_value));
    
    /*--- Create the right View for adding the matrices ---*/
    /* sigmaPoint is only pointing to one column
     * It already has incCol=0, so if we set nC=nState
     * it looks like a vector repeated nState times  */
    sigmaPoint->nC = nState;
    sigmaPoint->data = mean->data;
    
    /*--- Use two views for each half of SigmaPoints Matrix ---*/
    /* Use nState x nState Matrix: HelperMatrix */
    
    /* Set Sigma Points 1 to nState */
    /* SigmaPoints = x - sqrt(weight) P_sqrt*/
    HelperMatrix->data = SigmaPoints->data;
    result = lah_matAdd(HelperMatrix, 1, sigmaPoint, sqrt(ukal->weight), ukal->RootCovariance_state);
    
    /*Set Sigma Points nState to 2*nState */
    /* SigmaPoints = x - sqrt(weight) P_sqrt*/
    HelperMatrix->data += nState * SigmaPoints->incCol;
    result += lah_matAdd(HelperMatrix, 1, sigmaPoint, (-1.0) * sqrt(ukal->weight), ukal->RootCovariance_state);

    /* Set back sigmaPoin to old state*/
    sigmaPoint->nC = 1;
    sigmaPoint->data = SigmaPoints->data;
    return result;
}

