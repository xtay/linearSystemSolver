/****************************************************
some basic Matrix/vector operations that could be involved by all of the codes.

only symmetric matrix are supportted.

two kinds of features should be metioned here
    1.mpi support(next step)
    2.compressed storage
****************************************************/

#ifndef XTAY_LSS_UTIL_H
#define XTAY_LSS_UTIL_H

/**************************
  status tag definitions
  "LSS" stands for linear system solver
**************************/
#include <mpi.h>

#define LSS_FILE_DATA_LENGTH 10

#define LSS_SUCCESS 666

#define LSS_ERR_INPUT 1000
#define LSS_ERR_MEM 1001
#define LSS_ERR_ARG 1002
#define LSS_ERR_ITER 1003

//tags for matrix_vector_product() function
#define NEW_SPACE 1
#define OLD_SPACE 0

//something for parallel implementation
#define BLOCK_BASE(id, p, n) ((id)*(n)/(p))
#define BLOCK_SIZE(id, p, n) (BLOCK_BASE((id)+1, p, n) - BLOCK_BASE(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p)*(index+1)-1)/(n))

#define DATA_READ_TAG 2000
#define DATA_PRINT_TOKEN 3000
#define DATA_TRANSFER_TAG 4000
#define DATA_RETRIEVE_TAG 5000

#define USLEEP_TIME 0

/**************************
  some data structures
**************************/

typedef struct generalVector{
    double *array;
    int length;
}gVector;

typedef struct compressedMatrix{
    //a bunch of gVectors for each row of the matrix
    gVector *Lines;
    //size of the matrix;
    int nRows, nCols;
    int halfBandWidth;
}cMatrix;

typedef struct matrix_vector_product_env{
    int idMin;
    int idMax;
    int local_length;
    int remote_length;
    double *remote_array;
    double *remote_resArray;
    double *resBuffer;
    MPI_Request *sendRequest;
}mvpEnv;
/**************************
memory management
**************************/
int free_cMatrix(cMatrix *);
int free_gVector(gVector *);

/**************************
operations with element in the matrix 
**************************/
int attach_gVector_elem(gVector*, double*);
int remove_gVector_elem(gVector*, int);
int modify_gVector_elem(gVector*, int, double*);

int insert_cMatrix_elem(cMatrix*, int, int, double*);
int remove_cMatrix_elem(cMatrix*, int, int);
int modify_cMatrix_elem(cMatrix*, int, int, double*);
int get_cMatrix_elem(cMatrix*, int, int, double*);

/**************************
create an initial matrix/vector with data in a file
some memory space will be allocated within all these functions
**************************/
int fread_cMatrix(FILE*, cMatrix*);
int dup_cMatrix(cMatrix* dest, cMatrix* src);
int dup_gVector(gVector* dest, gVector* src);
int fread_gVector(FILE*, gVector*);
int pack_gVector(double *array, int length, gVector *res);

/**************************
print stuffs
**************************/
int print_cMatrix(cMatrix*);
int print_gVector(gVector*);

/**************************
basic calculation
**************************/
//a tag is used to tell the function to allocated space for the array or not
int set_mvp_env(mvpEnv *env, int halfBandWidth, int length);
int free_mvp_env(mvpEnv *);
int matrix_vector_product(cMatrix* cMat, gVector* gVec, gVector* res, int tag, mvpEnv *);
// a vector space will be allocated within this function
int inner_product(gVector* pVec1, gVector* pVec2, double* res);

#endif
