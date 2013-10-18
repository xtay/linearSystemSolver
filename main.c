#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "util/util.h"
#include "bicgstab/bicgstab.h"

int main(int argc, char *argv[])
{
    int lssStatus;
    int local_size;
    double innerProductRes;
    cMatrix cMat;
    gVector gVec;
    gVector xVec;
    gVector res_gVec;
    mvpEnv env;
    FILE *fMatPtr;
    FILE *fVecPtr;

    //mpi
    double elapsed_time;
    int id, nProc;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if(id == nProc - 1){
        if(argc != 3){
            printf("re-arange your arguments...\n");
            return LSS_ERR_ARG;
        }
        fMatPtr = fopen(argv[1], "r");
        fVecPtr = fopen(argv[2], "r");
    }

    if(id == 0) printf("test of read in compressed matrix\n");

    elapsed_time = -MPI_Wtime();
    lssStatus = fread_cMatrix(fMatPtr, &cMat);
    elapsed_time += MPI_Wtime();

    if(id == 0) printf("read time %10.6lf\n", elapsed_time);

    if(lssStatus == LSS_SUCCESS){

        elapsed_time = -MPI_Wtime();
        print_cMatrix(&cMat);
        elapsed_time += MPI_Wtime();

        if(id == 0) printf("print time %10.6lf\n", elapsed_time);
        if(id == 0) printf("SUCCESS: compressed matrix successfully read\n");
    }
    if(id == 0) printf("\n");

    if(id == 0) printf("test of read in Vector\n");
    lssStatus = fread_gVector(fVecPtr, &gVec);
    if(lssStatus == LSS_SUCCESS){
        print_gVector(&gVec);
        if(id == 0) printf("SUCCESS: vector successfully read\n\n");
    }

    if(id == 0) printf("test of duplicate a Vector\n");
    lssStatus = dup_gVector(&res_gVec, &gVec);
    if(lssStatus == LSS_SUCCESS){
        print_gVector(&res_gVec);
        if(id == 0) printf("SUCCESS: vector successfully duplicated\n\n");
    }

    if(id == 0) printf("test of inner_product of two Vectors\n");
    lssStatus = inner_product(&gVec, &gVec, &innerProductRes);
    if(lssStatus == LSS_SUCCESS){
        if(id == 0) printf("%10.3lf\n", innerProductRes);
        if(id == 0) printf("SUCCESS: vector inner_product success\n\n");
    }

    if(id == 0) printf("test of matrix_vector product\n");
    //res_gVec.array = (double *)calloc(gVec.length, sizeof(double));
    
    set_mvp_env(&env, cMat.halfBandWidth, gVec.length);
    lssStatus = matrix_vector_product(&cMat, &gVec, &res_gVec, NEW_SPACE, &env);
    free_mvp_env(&env);

    if(lssStatus == LSS_SUCCESS){
        print_gVector(&res_gVec);
        if(id == 0) printf("SUCCESS: matrix vector product success\n\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(id == 0) printf("test of solve a linear system by bicgstab method\n");
    lssStatus = bicgstab_solver(&cMat, &gVec, &xVec);
    if(lssStatus == LSS_SUCCESS){
        print_gVector(&xVec);
        if(id == 0) printf("SUCCESS: linear system successfully solved by bicgstab method\n\n");
    }
/*
*/

    MPI_Finalize();
    return 0;
}
