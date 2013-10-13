#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

/**************************
print stuffs
**************************/
int free_cMatrix(cMatrix* pMat)
{
    int m, n;
    int i, j, k;

    m = pMat->nRows;
    n = pMat->nCols;

    i = 0;
    while(i < m){
        free(pMat->Lines[i].array);
        i++;
    }
    free(pMat->Lines);

    return LSS_SUCCESS;
}

int free_gVector(gVector* pVec)
{
    free(pVec->array);
    return LSS_SUCCESS;
}

/**************************
operations with element in the matrix 
**************************/
int attach_gVector_elem(gVector* pVec, double* x);
int remove_gVector_elem(gVector* pVec, int index);
int modify_gVector_elem(gVector* pVec, int index, double* x);

int insert_cMatrix_elem(cMatrix*, int, int, double*);
int remove_cMatrix_elem(cMatrix*, int, int);
int modify_cMatrix_elem(cMatrix*, int, int, double*);
int get_cMatrix_elem(cMatrix*, int, int, double*);

/**************************
create an initial matrix/vector with data in a file
some memory space will be allocated within all these functions
**************************/

//assume that the file contains data in a specific form
//first line contains the size of the matrix
//rest of the file contains lines of nesscery numbers for each row of the matrix
int fread_cMatrix(FILE* fptr, cMatrix* pMat)
{
    int nRows, nCols; //tmp var for storing the size of the matrix;
    gVector *Lines; //tmp pointer points to the beginning of rows;
    
    char *tmpStr; //tmp string that stores a line of characters in the file
    char *token; //used for strtok()

    double *tmpData; //tmp pointer points to an array of numbers in a row;

    int i, j, k; //common indices

    fscanf(fptr, "%d", &nRows);
    fscanf(fptr, "%d", &nCols);
    pMat->nRows = nRows;
    pMat->nCols = nCols;

    //read in the size of the matrix
    //and allocate space for the Lines
    if(nRows && nCols){
        //an array that points to each row of the matrix
        Lines = (gVector *)malloc(sizeof(gVector) * nRows);
        tmpData = (double *)malloc(sizeof(double)*nRows);
        tmpStr = (char *)malloc(sizeof(char) * nRows * LSS_FILE_DATA_LENGTH);

        if(Lines == NULL || tmpData == NULL || tmpStr == NULL){
            printf("unable to allocate memory for Lines and temp spaces\n");
            return LSS_ERR_MEM;
        }

    }
    else{
        printf("Irregular matrix, exit!\n");
        return LSS_ERR_INPUT;
    }

    //read in the elements of the matrix
    if(Lines == NULL){
        printf("unable to allocate memory for tmpData\n");
        return LSS_ERR_MEM;
    }

    fgets(tmpStr, 10, fptr); //read the remaining \n of the first line

    i = 0;
    while(i < nRows){
        //obtaining data in a string
        fgets(tmpStr, sizeof(double) * nRows * LSS_FILE_DATA_LENGTH, fptr);

        //interpret the string into several double numbers
        j = 0; 
        token = strtok(tmpStr, " ");
        while(token != NULL){
            tmpData[j] = atof(token);
            token = strtok(NULL, " ");

            j++;
        }

        if(j > nCols - i){ //(nCols - i) is the number of elements that current compressed line could contain.
            printf("Irregular matrix, exit!\n");
            return LSS_ERR_INPUT;
        }

        //store to the right place
        Lines[i].array = malloc(sizeof(double)*j);
        Lines[i].length = j;
        memcpy(Lines[i].array, tmpData, sizeof(double)*j);

        i++;
    } 

    pMat->Lines = Lines;

    free(tmpData);
    free(tmpStr);
    return LSS_SUCCESS;
}

int dup_cMatrix(cMatrix* pDestMat, cMatrix* pSrcMat);

int dup_gVector(gVector* pDestVec, gVector* pSrcVec)
{
    int length;
    double *array;
    length = pSrcVec->length;

    array = (double *)malloc(sizeof(double)*length);
    if(array == NULL){
        printf("unable to allocate memory for array\n");
        return LSS_ERR_MEM;
    }

    memcpy(array, pSrcVec->array, sizeof(double)*length);

    pDestVec->array = array;
    pDestVec->length = length;
    return LSS_SUCCESS;
}

int fread_gVector(FILE* fptr, gVector* pVec)
{
    int length;
    double *array;

    int i, j, k;

    fscanf(fptr, "%d", &length);

    if(length == 0){
        printf("Irregular vector...\n");
        return LSS_ERR_INPUT;
    }

    array = (double *)malloc(sizeof(double)*length);

    i = 0;
    while(i < length){
        fscanf(fptr, "%lf", array+i);
        i++;
    }

    pVec->array = array;
    pVec->length = length;

    return LSS_SUCCESS;
}

int pack_gVector(double *array, int length, gVector *res)
{
    res->array = array;
    res->length = length;
    return LSS_SUCCESS;
}
/**************************
print stuffs
**************************/
//print it as a square matrix not a compressed one
//this is useless while process the solver alg.
//so, only a upper triangle symmetric part will be printed
int print_cMatrix(cMatrix* pMat)
{
    int nRows, nCols;

    double *tmpData;
    int length;
    
    int i, j, k;
    
    if(pMat == NULL){
        printf("Matrix seems not allocated\n");
        return LSS_ERR_MEM;
    }

    nRows = pMat->nRows;
    nCols = pMat->nCols;
    
    if(nRows == 0 || nCols == 0){
        printf("Cannot print a matrix with 0 rows or 0 cols!\n");
        return LSS_ERR_INPUT;
    }

    i = 0;
    while(i < nRows){
        tmpData = pMat->Lines[i].array;
        length = pMat->Lines[i].length;
        j = 0;
        while(j < i){
            printf("%10.3lf", 0.0);
            j++;
        }

        k = 0;
        while(k < length){
            printf("%10.3lf", tmpData[k]);
            k++;
        }

        j += k;
        while(j < nCols){
            printf("%10.3lf", 0.0);
            j++;
        }
        
        printf("\n");
        i++;
    }
    return LSS_SUCCESS;
}

int print_gVector(gVector* pVec)
{
    int length;
    double *array;
    
    int i, j, k;

    if(pVec == NULL){
        printf("Vector seems not allocated\n");
        return LSS_ERR_MEM;
    }

    length = pVec->length;
    array = pVec->array;

    if(length == 0){
        printf("Cannot print a vector with 0 elements!\n");
        return LSS_ERR_INPUT;
    }

    i = 0;
    while(i < length){
        printf("%10.3f", array[i]);
        i++;
    }
    printf("\n");
    
    return LSS_SUCCESS;
}

/**************************
basic calculation
**************************/
//note that this matrix/vector operations involves compressed matrix,
//so we should not call the inner_product() method below to deal with it.
int matrix_vector_product(cMatrix* pMat, gVector* pVec, gVector* res, int tag)
{
    double *resArray;
    double *arrayInMat;
    double *arrayInVec;

    int resLength;
    int lengthInMat;

    int nRows;
    int nCols;
    int length;

    int i, j, k;
    
    if(pMat == NULL || pVec == NULL){
        printf("Matrix or Vectors seems not allocated\n");
        return LSS_ERR_MEM;
    }

    length = pVec->length;
    nRows = pMat->nRows;
    nCols = pMat->nCols;

    if(length == 0 || nRows == 0 || nCols == 0){
        printf("Matrix or Vector has size 0\n");
        return LSS_ERR_INPUT;
    }
    if(length != nCols){
        printf("Matrix and Vector size mismatch\n");
        return LSS_ERR_INPUT;
    }

    resLength = nCols;
    if(tag)
        resArray = (double *)calloc(resLength, sizeof(double));
    else{
        resArray = res->array;
        //it's so important to reset the memory...
        memset(resArray, '\0', resLength*sizeof(double));
    }
    
    i = 0;
    arrayInVec = pVec->array;
    while(i < nCols){
        arrayInMat = pMat->Lines[i].array;
        lengthInMat = pMat->Lines[i].length;

        //elements on the diag
        resArray[i] += arrayInMat[0] * arrayInVec[i];
        j = 1;
        while(j < lengthInMat){
            //upper triangle part
            resArray[i] += arrayInMat[j] * arrayInVec[j+i];
            //lower triangle part
            resArray[j+i] += arrayInMat[j] * arrayInVec[i];
            j++;
        }

        i++;
    }
    res->length = resLength;
    res->array = resArray;

    return LSS_SUCCESS;
}

// a vector space will be allocated within this function
int inner_product(gVector* pVec1, gVector* pVec2, double* res)
{
    double result;
    double *array1, *array2;
    int length;
    int i, j, k;
    
    if(pVec1 == NULL || pVec2 == NULL){
        printf("Vectors seems not allocated\n");
        return LSS_ERR_MEM;
    }

    length = pVec1->length;
    if(length != pVec2->length){
        printf("two vectors of different size cannot do this...\n");
        return LSS_ERR_INPUT;
    }
    if(length == 0){
        printf("vector with 0 size cannot do this...\n");
        return LSS_ERR_INPUT;
    }
    array1 = pVec1->array;
    array2 = pVec2->array;

    i = 0;
    result = 0.0;
    while(i < length){
        result += (array1[i] * array2[i]);
        i++;
    }
    *res = result;

    return LSS_SUCCESS;
}
