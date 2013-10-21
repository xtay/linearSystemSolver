#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
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
    int halfBandWidth;
    gVector *Lines; //tmp pointer points to the beginning of rows;
    
    char *tmpStr; //tmp string that stores a line of characters in the file
    char *token; //used for strtok()

    double *tmpData; //tmp pointer points to an array of numbers in a row;

    int i, j, k, l; //common indices

    int id, nProc; //mpi
    int local_rows; //mpi
    MPI_Status status;
    

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    //read in the size of the matrix
    if(id == (nProc - 1)){
        fscanf(fptr, "%d", &nRows);
        fscanf(fptr, "%d", &nCols);
    }

    MPI_Bcast(&nRows, 1, MPI_INT, nProc-1, MPI_COMM_WORLD);
    if(!nRows)
        MPI_Abort(MPI_COMM_WORLD, LSS_ERR_INPUT);
    MPI_Bcast(&nCols, 1, MPI_INT, nProc-1, MPI_COMM_WORLD);

    pMat->nRows = nRows;
    pMat->nCols = nCols;

    //local size of the matrix
    local_rows = BLOCK_SIZE(id, nProc, nRows);

    //and allocate space for the Lines
    //an array that points to each row of the matrix
    Lines = (gVector *)malloc(sizeof(gVector) * local_rows);
    if(Lines == NULL){
        printf("unable to allocate memory for Lines in proc %d\n", id);
        return LSS_ERR_MEM;
    }

    //read in the elements of the matrix
    tmpData = (double *)malloc(sizeof(double)*nRows);
    if(tmpData == NULL){
        printf("unable to allocate memory for Lines and temp spaces, proc %d\n", id);
        return LSS_ERR_MEM;
    }

    if(id == (nProc - 1)){
        halfBandWidth = 0;
        
        tmpStr = (char *)malloc(sizeof(char) * nRows * LSS_FILE_DATA_LENGTH);

        if(tmpStr == NULL){
            printf("unable to allocate memory for Lines and temp spaces\n");
            return LSS_ERR_MEM;
        }


        fgets(tmpStr, 10, fptr); //read the remaining \n of the first line

        i = 0;
        while(i < nProc-1){

            //read and send data to some proccessor
            l = BLOCK_SIZE(i, nProc, nRows);
            k = 0;
            while(k < l){
                //obtaining data in a string
                fgets(tmpStr, sizeof(double) * nRows * LSS_FILE_DATA_LENGTH, fptr);

                //interpret the string into several double numbers
                j = 1; 
                token = strtok(tmpStr, " ");
                while(token != NULL){
                    tmpData[j] = atof(token);
                    token = strtok(NULL, " ");

                    j++;
                }
                tmpData[0] = j - 1; 

                if(j-1 > nCols - BLOCK_BASE(i, nProc, nRows) - k){ //(nCols - i) is the number of elements that current compressed line could contain.
                    printf("Irregular matrix, exit!, proc %d, row %d\n", id, k);
                    return LSS_ERR_INPUT;
                }

                halfBandWidth = halfBandWidth>(j-1) ? halfBandWidth : (j-1);

                //first element in tmpData tell the receiver how many data it will get..., the rest of tmpData is real data
                //argument tag is used so that the reciever would get the data in order
                MPI_Send(tmpData, j, MPI_DOUBLE, i, DATA_READ_TAG + k, MPI_COMM_WORLD);
                k++;
            }
            i++;
        } 

        //read data for this processor itself
        l =  local_rows;
        k = 0;
        while(k < l){
            fgets(tmpStr, sizeof(double) * nRows * LSS_FILE_DATA_LENGTH, fptr);

            //interpret the string into several double numbers
            j = 0; 
            token = strtok(tmpStr, " ");
            while(token != NULL){
                tmpData[j] = atof(token);
                token = strtok(NULL, " ");

                j++;
            }

            if(j > nCols - BLOCK_BASE(i, nProc, nRows) - k){ //(nCols - i) is the number of elements that current compressed line could contain.
                printf("Irregular matrix, exit!, proc %d, row %d\n", id, k);
                return LSS_ERR_INPUT;
            }
            
            halfBandWidth = halfBandWidth>(j-1) ? halfBandWidth : (j-1);
            //store to the right place
            Lines[k].length = j;
            Lines[k].array = malloc(sizeof(double)*j);
            memcpy(Lines[k].array, tmpData, sizeof(double)*j);

            k++;
        }
        pMat->Lines = Lines;

        free(tmpStr);
    }
    else{
        l = local_rows;
        k = 0;
        while(k < l){
            MPI_Recv(tmpData, nRows, MPI_DOUBLE, nProc-1, DATA_READ_TAG + k, MPI_COMM_WORLD, &status);
            j = tmpData[0];
            Lines[k].length = j;
            Lines[k].array = malloc(sizeof(double)*j);
            memcpy(Lines[k].array, tmpData+1, sizeof(double)*j);
            k++;
        }
        pMat->Lines = Lines;
    }

    MPI_Bcast(&halfBandWidth, 1, MPI_INT, nProc-1, MPI_COMM_WORLD);
    pMat->halfBandWidth = halfBandWidth;

    free(tmpData);
    return LSS_SUCCESS;
}

int dup_cMatrix(cMatrix* pDestMat, cMatrix* pSrcMat);

int dup_gVector(gVector* pDestVec, gVector* pSrcVec)
{
    int length;
    double *array;
    //mpi
    int nProc, id;
    int local_length;

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    length = pSrcVec->length;

    local_length = BLOCK_SIZE(id, nProc, length);

    array = (double *)malloc(sizeof(double)*local_length);

    if(array == NULL){
        printf("unable to allocate memory for array\n");
        return LSS_ERR_MEM;
    }

    memcpy(array, pSrcVec->array, sizeof(double)*local_length);

    pDestVec->array = array;
    pDestVec->length = length;
    return LSS_SUCCESS;
}

int fread_gVector(FILE* fptr, gVector* pVec)
{
    int length;
    double *array;
    double *tmpData;

    int i, j, k, l;

    //mpi
    int nProc, id;
    int local_length;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);

    if(id == nProc - 1)
        fscanf(fptr, "%d", &length);

    MPI_Bcast(&length, 1, MPI_INT, nProc-1, MPI_COMM_WORLD);

    if(length == 0){
        printf("Irregular vector...\n");
        MPI_Abort(MPI_COMM_WORLD, LSS_ERR_INPUT);
    }

    pVec->length = length;
    local_length = BLOCK_SIZE(id, nProc, length);

    array = (double *)malloc(sizeof(double)*local_length);

    //read for other processor
    if(id == nProc - 1){
        tmpData = (double *)malloc(sizeof(double)* (length/nProc+1));

        i = 0;
        while(i < nProc-1){
            l = BLOCK_SIZE(i, nProc, length);
            k = 0;
            while(k < l){
                fscanf(fptr, "%lf", tmpData+k);
                k++;
            }
            MPI_Send(tmpData, l, MPI_DOUBLE, i, DATA_READ_TAG, MPI_COMM_WORLD);
            i++;
        }
        l = local_length;
        k = 0;
        while(k < l){
            fscanf(fptr, "%lf", array+k);
            k++;
        }
        pVec->array = array;

        free(tmpData);
    }
    else{
        MPI_Recv(array, local_length, MPI_DOUBLE, nProc-1, DATA_READ_TAG, MPI_COMM_WORLD, &status);
        pVec->array = array;
    }

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

    //MPI
    int nProc, id;
    int local_rows;
    int print_token;//only processor who has the token could print
    MPI_Status status;
    
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

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    local_rows = BLOCK_SIZE(id, nProc, nRows);

    print_token = id;

    if(id > 0)
        MPI_Recv(&print_token, 1, MPI_INT, id-1, DATA_PRINT_TOKEN, MPI_COMM_WORLD, &status);

    print_token++;
    
    fflush(stdout);
    i = 0;
    while(i < local_rows){
        tmpData = pMat->Lines[i].array;
        length = pMat->Lines[i].length;
        j = 0;
        while(j < i+BLOCK_BASE(id, nProc, nRows)){
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
    //printf("%d\n", pMat->halfBandWidth);
    fflush(stdout);
    
    //in order to let the strings finish print on the screen, or the multi-process print will make the output a mess...
    //actually, the parameter is set to 0, but it works just fine...
    usleep(USLEEP_TIME);

    if(id < nProc-1)
        MPI_Send(&print_token, 1, MPI_INT, (id+1), DATA_PRINT_TOKEN, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    usleep(USLEEP_TIME);

    return LSS_SUCCESS;
}

int print_gVector(gVector* pVec)
{
    int length;
    double *array;
    
    int i, j, k;

    //mpi
    int nProc, id;
    int local_length;
    int print_token;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

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

    local_length = BLOCK_SIZE(id, nProc, length);

    if(id > 0)
        MPI_Recv(&print_token, 1, MPI_INT, id-1, DATA_PRINT_TOKEN, MPI_COMM_WORLD, &status);

    print_token++;
    fflush(stdout);
    i = 0;
    while(i < local_length){
        printf("%10.3f", array[i]);
        i++;
    }
    fflush(stdout);
    //in order to let the strings finish print on the screen, or the multi-process print will make the output a mess...
    usleep(USLEEP_TIME);

    if(id < nProc-1)
        MPI_Send(&print_token, 1, MPI_INT, (id+1), DATA_PRINT_TOKEN, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    usleep(USLEEP_TIME);

    if(id == 0)
        printf("\n");
    
    return LSS_SUCCESS;
}

/**************************
basic calculation
**************************/
//note that this matrix/vector operations involves compressed matrix,
//so we should not call the inner_product() method below to deal with it.
int set_mvp_env(mvpEnv *env, int halfBandWidth, int length){
    int local_length;
    int remote_length;
    int idMax, idMin, idWidth;
    int nProc, id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);

    idMin = BLOCK_OWNER(BLOCK_BASE(id, nProc, length) - halfBandWidth + 1, nProc, length);
    idMax = BLOCK_OWNER(BLOCK_BASE(id+1, nProc, length)-1 + halfBandWidth - 1, nProc, length);
    
    local_length = BLOCK_SIZE(id, nProc, length); 
    remote_length = BLOCK_BASE(idMax + 1, nProc, length) - BLOCK_BASE(id+1, nProc, length);

    idWidth = id - idMin > idMax - id ? id - idMin : idMax - id;

    env->idMin = idMin;
    env->idMax = idMax;
    env->local_length = local_length;
    env->remote_length = remote_length;
    env->remote_array = (double *)calloc(remote_length, sizeof(double));
    env->remote_resArray = (double *)calloc(remote_length, sizeof(double));
    env->resBuffer = (double *)calloc(local_length, sizeof(double));
    env->sendRequest = (MPI_Request *)calloc(idWidth, sizeof(MPI_Request));

    return LSS_SUCCESS;
}

int free_mvp_env(mvpEnv *env){
    free(env->remote_array);
    free(env->remote_resArray);
    free(env->resBuffer);
    free(env->sendRequest);
    return LSS_SUCCESS;
}
int matrix_vector_product(cMatrix* pMat, gVector* pVec, gVector* res, int tag, mvpEnv *env)
{
    double *arrayInMat;
    double *arrayInVec;

    int remote_length;
    int lengthInMat;
    int lengthInVec;
    int halfBandWidth;

    int nRows;
    int nCols;
    int length;

    int i, j, k;

    //mpi
    double *remote_array;
    double *local_resArray;
    double *remote_resArray;
    double *resBuffer;
    int idMin, idMax, idWidth;
    int nProc, id;
    int local_rows;
    int local_length;
    MPI_Status status;
    MPI_Request *sendRequest;
    int target;

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

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

    local_rows = BLOCK_SIZE(id, nProc, nRows);

    arrayInVec = pVec->array;
    lengthInVec = pVec->length;

    halfBandWidth = pMat->halfBandWidth;

    idMin = env->idMin;
    idMax = env->idMax;
    local_length = env->local_length;
    remote_length = env->remote_length;
    remote_array = env->remote_array;
    remote_resArray = env->remote_resArray;
    resBuffer = env->resBuffer;
    sendRequest = env->sendRequest;

    memset(remote_resArray, '\0', remote_length*sizeof(double));

    if(tag)
        local_resArray = (double *)calloc(local_length, sizeof(double));
    else{
        local_resArray = res->array;
        memset(local_resArray, '\0', local_length*sizeof(double));
    }
    //retrieve part of the vector which is nesscery to perform a local calculation
    while(idMin < id){
        target = idMin<0 ? MPI_PROC_NULL : idMin;
        MPI_Isend(arrayInVec, local_length, MPI_DOUBLE, target, DATA_TRANSFER_TAG, MPI_COMM_WORLD, sendRequest+id-idMin);
        idMin++;
    }

    //finish the local calculation
    i = 0;
    while(i < local_length){
        arrayInMat = pMat->Lines[i].array;
        lengthInMat = pMat->Lines[i].length;

        //elements on the diag
        local_resArray[i] += arrayInMat[0] * arrayInVec[i];
        j = 1;
        while(i + j < local_length && j < lengthInMat){
            //upper triangle part
            local_resArray[i] += arrayInMat[j] * arrayInVec[j+i];
            //lower triangle part
            local_resArray[j+i] += arrayInMat[j] * arrayInVec[i];
            j++;
        }

        i++;
    }

    //
    while(idMax > id){
        target = idMax>(nProc-1) ? MPI_PROC_NULL : idMax;
        MPI_Recv(remote_array + BLOCK_BASE(idMax, nProc, length) - BLOCK_BASE(id+1, nProc, length), BLOCK_SIZE(idMax, nProc, length), MPI_DOUBLE, target, DATA_TRANSFER_TAG, MPI_COMM_WORLD, &status);
        idMax--;
    }


    //do the rest of the calculation
    i = 0;
    while(i < local_length){
        arrayInMat = pMat->Lines[i].array;
        lengthInMat = pMat->Lines[i].length;

        j = local_length - i;
        while(j < lengthInMat){
            //upper triangle part
            local_resArray[i] += arrayInMat[j] * remote_array[j+i - local_length];
            //lower triangle part
            remote_resArray[j+i - local_length] += arrayInMat[j] * arrayInVec[i];
            j++;
        }
        i++;
    }


    //exchange the answers
    idMin = BLOCK_OWNER(BLOCK_BASE(id, nProc, length) - halfBandWidth + 1, nProc, length);
    idMax = BLOCK_OWNER(BLOCK_BASE(id+1, nProc, length)-1 + halfBandWidth - 1, nProc, length);

    while(idMax > id){
        target = idMax>(nProc-1) ? MPI_PROC_NULL : idMax;
        MPI_Isend(remote_resArray + BLOCK_BASE(idMax, nProc, length) - BLOCK_BASE(id+1, nProc, length), BLOCK_SIZE(idMax, nProc, length), MPI_DOUBLE, target, DATA_TRANSFER_TAG, MPI_COMM_WORLD, sendRequest+idMax-id);
        idMax--;
    }

    while(idMin < id){
        target = idMin<0 ? MPI_PROC_NULL : idMin;
        MPI_Recv(resBuffer, local_length, MPI_DOUBLE, target, DATA_TRANSFER_TAG, MPI_COMM_WORLD, &status);
        i = 0;
        while(i < local_length){
            local_resArray[i] += resBuffer[i];
            //printf("id%d, i%d, val%lf\n", id, i, resArray[i]);
            i++;
        }
        idMin++;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    res->length = nCols;
    res->array = local_resArray;

    return LSS_SUCCESS;
}

// a vector space will be allocated within this function
int inner_product(gVector* pVec1, gVector* pVec2, double* res)
{
    double local_result;
    double *array1, *array2;
    int length;
    int i, j, k;

    //mpi
    int nProc, id;
    int local_length;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
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

    local_length = BLOCK_SIZE(id, nProc, length);

    array1 = pVec1->array;
    array2 = pVec2->array;

    i = 0;
    local_result = 0.0;
    while(i < local_length){
        local_result += (array1[i] * array2[i]);
        i++;
    }

    //
    if(id == 0){
        *res = local_result;
        i = 1;
        while(i < nProc){
            MPI_Recv(&local_result, 1, MPI_DOUBLE, i, DATA_TRANSFER_TAG, MPI_COMM_WORLD, &status);
            *res += local_result;
            i++;
        }
    }
    else{
        MPI_Send(&local_result, 1, MPI_DOUBLE, 0, DATA_TRANSFER_TAG, MPI_COMM_WORLD);
    }

    MPI_Bcast(res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return LSS_SUCCESS;
}
