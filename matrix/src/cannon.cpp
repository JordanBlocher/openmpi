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

#include "matlib.hpp"

#define ATAG 1
#define BTAG 2
#define CTAG 3
#define DIETAG 4
#define MAXNUMPROCS 20

//MPI Datatypes
MPI_Datatype MPI_Vector;

template <typename T>
void work(T*,T*,T*);

int N, numprocs, size, size2;
MPI_Request request;
MPI_Status status;
int rank, err=0;

int main (int argc, char *argv[]) {
    if(argc<2)
    {
        printf("Usage: ./bin/matrix <N> \n");
        return 0;
    }

    int i,j;
    N = atoi(argv[1]);
  
    std::chrono::high_resolution_clock::time_point startclock, stopclock;
    std::chrono::duration<double> time_span;
 
    // Init MPI & get numprocs and rank
    MPI_Init(&argc,&argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    CHECKMPI(err); 
    
    //if(N % numprocs != 0 || numprocs < 4 || fmod(sqrt(numprocs), 1) > 10e-6)
      //  return 0;

    size = N / (int)sqrt(numprocs);
    size2 = size*size;
    Matrix<int> *Ap, *Bp, *Cp;
    int *Arow, *Brow, *Crow;
    Arow = new int[size2];
    Brow = new int[size2];
    Crow = new int[size2];

    if(rank == 0)
    {
        Ap = new Matrix<int>(N,N);
        Bp = new Matrix<int>(N,N);
        Cp = new Matrix<int>(N,N);
        Cp->clear();

        Matrix<int> &A = *Ap;
        Matrix<int> &B = *Bp;

        // Initialize matrix A random, B identity
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                A(i, j) = i*N+j;
                if(i==j)
                    B(N-i-1, j) = 1;
                else B(N-i-1, j) = 0;
            }
        }
        A.setRowIndexing(size, N, true);
        B.setRowIndexing(size, N, true);

        pmatprintf(A, size, "A");
        pmatprintf(B, size, "B");
        
        startclock = std::chrono::high_resolution_clock::now();

        A.shift(size, N, false, true);
        B.shift(size, N, true, true);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0; i<N/size; ++i)
    {
        // Scatter B to processes
        err = MPI_Scatter(Bp->cmat, size2, MPI_INT, &(B.cmat[rank*size2]), size2, MPI_INT, 0, MPI_COMM_WORLD);
        CHECKMPI(err); 
        
        // Scatter A to processes
        err = MPI_Scatter(Bp->cmat, size2, MPI_INT, &(A.cmat[rank*size2]), size2, MPI_INT, 0, MPI_COMM_WORLD);
        CHECKMPI(err); 

        // Perform work on matrix
        work(&(A.cmat[rank*size2]), &(B.cmat[rank*size2]), &(C.cmat[rank*size2]));

        // Gather back result to C
        err= MPI_Gather(&(C.cmat[rank*size2]), size2, MPI_INT, Cp->cmat, size2, MPI_INT, 0, MPI_COMM_WORLD);
        CHECKMPI(err); 
        
/*
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0)
        {
        pmatprintf(A, size, "A");
        cmatprintf(A);
        pmatprintf(B, size, "B");
        cmatprintf(B);
        printf("Indexing ..\n");
        pmatprintf(C, size, "C");
        cmatprintf(C);
        }
*/
        
        if (i!=(N/size-1))
        {
            A.shift(size, N);
            B.shift(size, N, true);
        }
    }

    if(rank == 0)
    {
        stopclock = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stopclock - startclock);

        std::cout<<" = C\n";
        pmatprintf(C, size, "C");

        printf("\nProgram exec time: %f seconds \n\n",(double)time_span.count());
    }

    MPI_Finalize();
    return 0;
}

template <typename T>
void work(T* A, T* B, T* C)
{
   // submatprintf(A,size,size,"Awork", rank); std::cout<<"\n x \n "; submatprintf(B,size,size,"Bwork", rank);
//    std::cout<<"\n";
    int i,j,k;
    for(i=0;i<size;i++)
    {
   //     printf("row %d\n", i);
        for(j=0;j<size;j++)
        {
            for(k=0;k<size;k++) 
            {
                C[j+i*size]+=A[k+i*size]*B[j+k*size];
 //               printf("\tC(%d,%d) += A(%d,%d)*B(%d,%d)\n",i,j,i,k,k,j);
  //              printf("\t\t%d += %d*%d\n",C[j+i*size],A[k+i*size],B[j+k*size]);
            }
        }
    }
  //  std::cout<<" = \n"; submatprintf(C,size,size,"Cwork", rank);
  //  std::cout<<"\n\n";
}


