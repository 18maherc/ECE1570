//Code by Joey Black and Curtis Maher

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define CACHE_LINE_SIZE 128

void matrixMultiply(const int, const double *, const double *, double *);
void blockedMatrixMultiply(const int M, const double *A, const double *B, double *C);
void subMatrixMultiply(const int M, const double *A, const double *B, double *C,
                       const int rowStartA, const int columnStartB, const int columnStartA, 
                       const int nRowStartA, const int nColumnStartB, const int nColumnStartA);

int main()
{
    //Use clock for stats on algo and time for randomization
    clock_t start, end;
    time_t t;
    double cpu_time_used;

    const int size = 2000;
    double *A = (double *) malloc(size * size * sizeof(double));
    double *B = (double *) malloc(size * size * sizeof(double));
    double *C1 = (double *) malloc(size * size * sizeof(double));
    double *C2 = (double *) malloc(size * size * sizeof(double));

    //Initialize all the matrices
    srand((unsigned) time(&t));
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            A[size * i + j] = rand() % 1000;
            B[size * i + j] = rand() % 1000;
            C1[size * i + j] = 0;
            C2[size * i + j] = 0;
        }
    }

    //Run naive algorithm and get the time
    start = clock();
    matrixMultiply(size, A, B, C1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);
    
    //Run the block algorithm and get the time
    start = clock();
    blockedMatrixMultiply(size, A, B, C2);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);

    //Check for errors
    for(int i = 0; i < size * size; i++)
    {
        if(C1[i] != C2[i])
        {
            printf("fail\n");
            break;
        }
    }

    //Free all the memory
    free(A);
    free(B);
    free(C1);
    free(C2);
}

void matrixMultiply(const int M, const double *A, const double *B, double *C)
{
int i , j , k;
    for ( i = 0; i < M; ++i) {
        for ( j = 0; j < M; ++j) {
            for (k = 0; k < M; ++k)
                C[i*M + j] += A[i*M+k] * B[k*M+j];
        }
    }
}

void blockedMatrixMultiply(const int M, const double *A, const double *B, double *C)
{
    int i, j, k;
    int outsideValue = M % CACHE_LINE_SIZE;
    int cacheIndex = CACHE_LINE_SIZE;
    int rowASize, colBSize, colASize;
    
    //Make sure that the size is greater than the cache line size
    //If not, just use single block since it won't benefit from blocking
    if(M < cacheIndex)
    {
        cacheIndex = M;
        outsideValue = 0;
    }

    /*
        Goes through each block and does matrix multiplication
        A block is size CacheBlock x CacheBlock
        Handles matrices that aren't mutliples of block size
        This algo hovers around 15-20 times faster than Naive, something else is probably affected it
            but Im still willing to say its at least 10x faster
    */

    for(i = 0; i < M; i+= CACHE_LINE_SIZE)
    {
        for(k = 0; k < M; k+= CACHE_LINE_SIZE)
        {
            for(j = 0; j < M; j+= CACHE_LINE_SIZE)
            {
                rowASize = cacheIndex, colBSize = cacheIndex, colASize = cacheIndex;
                if(outsideValue != 0)
                {
                    if(i + CACHE_LINE_SIZE >= M && outsideValue != 0)
                        rowASize = outsideValue;
                    if(j + CACHE_LINE_SIZE >= M && outsideValue != 0)
                        colBSize = outsideValue;
                    if(k + CACHE_LINE_SIZE >= M && outsideValue != 0)
                        colASize = outsideValue;
                } 
                subMatrixMultiply(M, A, B, C, i, j, k, rowASize, colBSize, colASize);
            }
        }
    }
}

void subMatrixMultiply(const int M, const double *A, const double *B, double *C,
                       const int rowStartA, const int columnStartB, const int columnStartA, 
                       const int nRowStartA, const int nColumnStartB, const int nColumnStartA)
{
    //Goes through a block and does matrix multiplication
    //Goes row by row through A and B and changes the write address everytime
    //This should make cache access on B easier since its going by row and not
    //By column
    int i, j, k;
    for(i = rowStartA; i < (nRowStartA + rowStartA); i++)
    {
        for(k = columnStartA; k < (nColumnStartA + columnStartA); k++)
        {
            for(j = columnStartB; j < (nColumnStartB + columnStartB); j++)
            {
                C[i*M + j] += A[i*M + k] * B[k*M + j];
            }
        }
    }
}