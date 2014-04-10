#define CHECKMPI(err) if(err){ fprintf(stderr,"Communicator failed.\n");MPI_Abort(MPI_COMM_WORLD,1);};

#include <stdio.h>
#include <mpi.h>
#include "string.h"
#include <time.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <sstream>
#include <fstream>
#include <cmath>
#include <array>
#include <chrono>
#include <algorithm>

#define ATAG 1
#define BTAG 2
#define DIETAG 3
#define MAXNUMPROCS 20

template <typename T>
void work(T* A, T* B, T* C, int, int);
template <typename T>
void matprintf(T* A);
template <typename T>
void submatprintf(T* A, int, int, int, int);

int rows, cols, numprocs, size;
MPI_Request request;
MPI_Status status;
int rank, err=0;

inline int index(int i) {return (rows / numprocs)*(i);}

int main (int argc, char *argv[]) {
    if(argc<3)
    {
        printf("Usage: ./bin/matrix <rows> <cols> \n");
        return 0;
    }

    int i,j;
    rows = atoi(argv[1]);
    cols = atoi(argv[2]);
  
    std::chrono::high_resolution_clock::time_point startclock, stopclock;
    std::chrono::duration<double> time_span;
 
    // Init MPI & get numprocs and rank
    MPI_Init(&argc,&argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    CHECKMPI(err); 

    size = rows / numprocs;
    int extra = rows % numprocs;
    int* B = new int[rows*cols];
    int* A = new int[rows*cols];
    int* C = new int[rows*cols];

    if(rank == 0)
    {
       // Initialize matrix A random, B identity
        for(i=0; i<rows; i++)
        {
            for(j=0; j<cols; j++)
            {
                A[j + i * cols] = rand() % 20 + 1;
                if(i==j)
                    B[j + i * cols]= 1;
                else
                    B[j + i * cols]= 0;
                C[j + i * cols]= 0;
            }
        }
        std::cout<<"A\n";
        matprintf(A);
        std::cout<<"x B\n";
        matprintf(B);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        startclock = std::chrono::high_resolution_clock::now();

    // Broadcast B to processes
    err = MPI_Bcast(B, rows*cols, MPI_INT, 0, MPI_COMM_WORLD);
    CHECKMPI(err); 
    
    // Scatter A to processes
    err = MPI_Scatter(A, size*cols, MPI_INT, &A[index(rank)*cols], size*cols, MPI_INT, 0, MPI_COMM_WORLD);
    CHECKMPI(err); 

    // Perform work on matrix
    work(A, B, C, index(rank), index(rank+1));

    // Work on any extra rows
    if(rank == 0 && extra > 0)
    {
        work(A, B, C, rows-extra, rows);
    }
 
    // Gather back result to C
    err= MPI_Gather(&C[index(rank)*cols], size*cols, MPI_INT, C, size*cols, MPI_INT, 0, MPI_COMM_WORLD);
    CHECKMPI(err); 

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        stopclock = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stopclock - startclock);
        std::cout<<"= C\n";
        matprintf(C);
        printf("\nProgram exec time: %f seconds \n\n",(double)time_span.count());
    }

    MPI_Finalize();
    return 0;
}

template <typename T>
void work(T* A, T* B, T* C, int rs, int re)
{
    int i,j,k;
    for(i=rs;i<re;i++)
        for(j=0;j<cols;j++)
        {
           C[j + i * cols]= 0;
            for(k=0;k<cols;k++) 
            {
                C[j + i * cols]+=A[k + i * cols]*B[j + k * cols];
            }
        }
}

template <typename T>
void matprintf(T* A)
{
	int i,j,m,n;
    m = rows > 15 ? 15 : rows;
    n = cols > 15 ? 15 : cols;
	for(i=0;i<m;i++)
    {
        printf("\t ");
        for(j=0;j<n;j++)    
            printf("%d ", A[i * cols + j]);
        printf("\n");
    }
    printf("\n\n");
}
