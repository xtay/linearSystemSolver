#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "../util/util.h"

#define PRECISION 0.001
#define ITERATION_LIMIT 4848

//to know what are these variables stands for, see wiki page for "bicgstab"
int bicgstab_solver(cMatrix *pA, gVector *pb, gVector *px)
{
    gVector rVec;
    gVector rcVec;
    gVector p_Vec;
    gVector sVec;
    gVector vVec;
    gVector tVec;

    gVector tmpVec;

    //space for array b suppose to have been allocated
    double *bArray;

    //spaces for arrays below need to be allocated in this function
    double *xArray;

    double *rArray;
    double *rcArray;
    double *p_Array;
    double *sArray;
    double *vArray;
    double *tArray;

    double *tmpArray;

    int nRows;
    int nCols;
    int bLength;
    int xLength;
    int length;

    int i, j, k;
    double rho=1;
    double alpha=1;
    double beta;
    double omega=1;
    
    double rhoNew;
    double omegaNew;
    
    double tmpDouble1, tmpDouble2;

    double residualNorm;

    //mpi
    int nProc, id;
    int local_length;
    mvpEnv env;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);

    nRows = pA->nRows;
    nCols = pA->nCols;
    bLength = pb->length;
    xLength = nCols;
    length = xLength;
    local_length = BLOCK_SIZE(id, nProc, length);

    if(nRows == 0 || nCols == 0 || bLength == 0){
        printf("Failed: cannot solve an linear system with 0 size...\n");
        return LSS_ERR_INPUT;
    }

    if(nRows != bLength || xLength != bLength){
        printf("Failed: This linear system is trival...\n");
        return LSS_ERR_INPUT;
    }

    bArray = pb->array;
    //only the space for array x should be kept after this function finished
    xArray = (double *)calloc(local_length, sizeof(double));

    //spaces allocated for these arrays needs to be freed
    rArray = (double *)calloc(local_length, sizeof(double));
    vArray = (double *)calloc(local_length, sizeof(double));
        //these two vectors' space will be allocated by duplication later
        //rcArray = (double *)calloc(local_length, sizeof(double));
        //p_Array = (double *)calloc(local_length, sizeof(double));
    sArray = (double *)calloc(local_length, sizeof(double));
    tArray = (double *)calloc(local_length, sizeof(double));
    tmpArray = (double *)calloc(local_length, sizeof(double));

    //package these arrays into relevant gVectors
    pack_gVector(xArray, length, px);

    pack_gVector(rArray, length, &rVec);
    pack_gVector(vArray, length, &vVec);
    pack_gVector(sArray, length, &sVec);
    pack_gVector(tArray, length, &tVec);
    pack_gVector(tmpArray, length, &tmpVec);

    set_mvp_env(&env, pA->halfBandWidth, length);
    matrix_vector_product(pA, px, &tmpVec, OLD_SPACE, &env);
    tmpArray = tmpVec.array;

    i = 0;
    while(i < local_length){
        rArray[i] = bArray[i] - tmpArray[i];
        i++;
    }

    //here in this function a new space allocated for array in rcVec
    dup_gVector(&rcVec, &rVec);
    rcArray = rcVec.array;
    //here in this function a new space allocated for array in p_Vec
    dup_gVector(&p_Vec, &vVec);
    p_Array = p_Vec.array;

    i = 0;
    inner_product(&rVec, &rVec, &residualNorm);
    residualNorm = sqrt(residualNorm);

    while(residualNorm > PRECISION && i < ITERATION_LIMIT){
//        if(id == 1) printf("iteration: %d\n", i);

        inner_product(&rcVec, &rVec, &rhoNew);
        beta = (rhoNew/rho)*(alpha/omega);

        k = 0;
        while(k < local_length){
            p_Array[k] = rArray[k] + beta * (p_Array[k] - omega * vArray[k]);
            k++;
        }
//        if(id == 0) printf("p:\n");
//        print_gVector(&p_Vec);

        matrix_vector_product(pA, &p_Vec, &vVec, OLD_SPACE, &env);

        inner_product(&rcVec, &vVec, &tmpDouble1);
        alpha = rhoNew/tmpDouble1;

        k = 0;
        while(k < local_length){
            sArray[k] = rArray[k] - alpha * vArray[k];
            k++;
        }

        matrix_vector_product(pA, &sVec, &tVec, OLD_SPACE, &env);

        inner_product(&tVec, &sVec, &tmpDouble1);
        inner_product(&tVec, &tVec, &tmpDouble2);
        omegaNew = tmpDouble1 / tmpDouble2;

        k = 0;
        while(k < local_length){
            xArray[k] = xArray[k] + alpha * p_Array[k] + omegaNew * sArray[k];
            rArray[k] = sArray[k] - omegaNew * tArray[k];
            k++;
        }

        inner_product(&rVec, &rVec, &residualNorm);
        residualNorm = sqrt(residualNorm);

        omega = omegaNew;
        rho = rhoNew;

        i++;
    }
    
    free_mvp_env(&env);
    free(rArray);
    free(vArray);
    free(rcVec.array);
    free(p_Vec.array);
    free(sArray);
    free(tArray);
    free(tmpArray);

    if(i == 100){
        if(id == 0) printf("Failed: cannot get it solved within ITERATION LIMITATION\n");
        if(id == 0) printf("iteration times: %d\n residualNorm: %lf\n", i, residualNorm);
        return LSS_ERR_ITER;
    }
    else{
        if(id == 0)
            if(id == 0) printf("iteration times: %d\n residualNorm: %lf\n", i, residualNorm);
    }

    return LSS_SUCCESS;
}
