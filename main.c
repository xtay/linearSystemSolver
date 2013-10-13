#include <stdio.h>
#include <stdlib.h>

#include "util/util.h"
#include "bicgstab/bicgstab.h"

int main(int argc, char *argv[])
{
    int lssStatus;
    double innerProductRes;
    cMatrix cMat;
    gVector gVec;
    gVector xVec;
    gVector res_gVec;
    FILE *fMatPtr;
    FILE *fVecPtr;

    if(argc != 3){
        printf("re-arange your arguments...\n");
        return LSS_ERR_ARG;
    }
    fMatPtr = fopen(argv[1], "r");
    fVecPtr = fopen(argv[2], "r");

    printf("test of read in compressed matrix\n");
    lssStatus = fread_cMatrix(fMatPtr, &cMat);
    if(lssStatus == LSS_SUCCESS){
        print_cMatrix(&cMat);
        printf("SUCCESS: compressed matrix successfully read\n");
    }

    printf("\n");

    printf("test of read in Vector\n");
    lssStatus = fread_gVector(fVecPtr, &gVec);
    if(lssStatus == LSS_SUCCESS){
        print_gVector(&gVec);
        printf("SUCCESS: vector successfully read\n\n");
    }

    printf("test of duplicate a Vector\n");
    lssStatus = dup_gVector(&res_gVec, &gVec);
    if(lssStatus == LSS_SUCCESS){
        print_gVector(&res_gVec);
        printf("SUCCESS: vector successfully duplicated\n\n");
    }

    printf("test of inner_product of two Vectors\n");
    lssStatus = inner_product(&gVec, &gVec, &innerProductRes);
    if(lssStatus == LSS_SUCCESS){
        printf("%10.3lf\n", innerProductRes);
        printf("SUCCESS: vector inner_product success\n\n");
    }

    printf("test of matrix_vector product\n");
    res_gVec.array = (double *)calloc(gVec.length, sizeof(double));
    lssStatus = matrix_vector_product(&cMat, &gVec, &res_gVec, OLD_SPACE);
    if(lssStatus == LSS_SUCCESS){
        print_gVector(&res_gVec);
        printf("SUCCESS: matrix vector product success\n\n");
    }

    printf("test of solve a linear system by bicgstab method\n");
    lssStatus = bicgstab_solver(&cMat, &gVec, &xVec);
    if(lssStatus == LSS_SUCCESS){
        print_gVector(&xVec);
        printf("SUCCESS: linear system successfully solved by bicgstab method\n\n");
    }

    return 0;
}
