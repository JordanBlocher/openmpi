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
void master(Matrix<T>&, Matrix<T>&, Matrix<T>&);
template <typename T>
void slave(Matrix<T>&, Matrix<T>&, Matrix<T>&);
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
    Matrix<int> B(N, N);
    Matrix<int> A(N, N);
    Matrix<int> C(N, N);
    C.clear();
 
    if(rank == 0)
    {
       // Initialize matrix A random, B identity
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                A(i, j) = rand() % 20 + 1;
                C(i, j) = 0;
                if(i==j)
                    B(i, j) = 1;
                else B(i, j) = 0;
                printf("\t(%d,%d) ", i, j);
            }
            printf("\n");
        }
        B(0, 1) = 2;
        B(0, 2) = 3;
        std::cout<<"A\n";
        matprintf(A);
        std::cout<<"B\n";
        matprintf(B);
        std::cout<<"\n\n\n";
        //k=0;
        /*for(i=0; i<=N; i+=size)
            for(j=0; j<N; j+=size)
            {
                A.submat(i, j, size, k);
                submatprintf(A.smat[k], size, size, "Amaster");
                B.submat(i, j, size, k);
                submatprintf(B.smat[k], size, size, "Bmaster");
                C.submat(i, j, size, k);
                k++;
            }
*/
        startclock = std::chrono::high_resolution_clock::now();
       
    }
    A.partition(size, N, rank);
    B.partition(size, N, rank, true);
    if(rank == 0)
    {
    printf("A partitioned\n");
    matprintf(A);
    printf("B partitioned\n");
    matprintf(B);
    }
    C.partition(size, N, rank);

/*         
    if(rank == 0)
        master<int>(A, B, C);
    else
        slave<int>(A, B, C);
  */   
    printf("starting sends\n\n\n\n");
    // Scatter B to processes
    err = MPI_Scatter(B.mat, size2, MPI_INT, &(B.mat[rank*N]), size2, MPI_INT, 0, MPI_COMM_WORLD);
    CHECKMPI(err); 
    
    // Scatter A to processes
    err = MPI_Scatter(A.mat, size2, MPI_INT, &(A.mat[rank*N]), size2, MPI_INT, 0, MPI_COMM_WORLD);
    CHECKMPI(err); 

    // Perform work on matrix
    printf("rank %d doing work\n", rank);
    work(&(A.mat[rank*N]), &(B.mat[rank*N]), &(C.mat[rank*N]));

    // Gather back result to C
    err= MPI_Gather(&(C.mat[rank*N]), size2, MPI_INT, C.mat, size2, MPI_INT, 0, MPI_COMM_WORLD);
    CHECKMPI(err); 

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout<<"C done\n";
    matprintf(C);


    if(rank == 0)
    {
        stopclock = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stopclock - startclock);
        printf("\nProgram exec time: %f seconds \n\n",(double)time_span.count());
    }

    MPI_Finalize();
    return 0;
}

template <typename T>
void work(T* A, T* B, T* C)
{
    //submatprintf(A,size,size,"Awork", rank); std::cout<<" x \n "; submatprintf(B,size,size,"Bwork");
    int i,j,k;
    for(i=0;i<size;i++)
    {
     //   printf("row %d\n", i);
        for(j=0;j<size;j++)
        {
            C[j+i*size] = 0;
            for(k=0;k<size;k++) 
            {
                C[j+i*size]+=A[k+i*size]*B[j+k*size];
       //         printf("\tC(%d,%d) += A(%d,%d)*B(%d,%d)\n",i,j,i,k,k,j);
         //       printf("\t\t%d += %d*%d\n",C[j+i*size],A[k+i*size],B[j+k*size]);
            }
        }
    }
    //std::cout<<" = "; submatprintf(C,size,size,"Cwork");
    //std::cout<<"\n\n";
}

template <typename T>
void master(Matrix<T>& A, Matrix<T>& B, Matrix<T>& C)
{
    int i,j,k=1;
    for(i=0; i<size; ++i)
    {
        for(j=0; j<size; ++j)
        {
            if(i==0 && j==0) continue;
            printf("sending\n");
            err = MPI_Send(A.smat[k], size2, MPI_INT, k, ATAG, MPI_COMM_WORLD);
            CHECKMPI(err); 
            err = MPI_Send(B.smat[k], size2, MPI_INT, k, BTAG, MPI_COMM_WORLD);
            CHECKMPI(err); 
            k++;
        }
    }
    
    work(A.smat[0], B.smat[0], C.smat[0]);
    k=1;
    for(i=0; i<size; ++i)
    {
        for(j=0; j<size; ++i)
        {
            if(i==0 && j==0) continue;
            std::cout<<"waiting...\n";
            err = MPI_Recv(C.smat[k], size2, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            std::cout<<"recieved\n";
            CHECKMPI(err); 
            err = MPI_Send(0, 0, MPI_INT, k, DIETAG, MPI_COMM_WORLD);
            CHECKMPI(err); 
            k++;
        }
    }
    std::cout<<"C done\n";
    matprintf(C);

    return;
}

template <typename T>
void slave(Matrix<T>& A, Matrix<T>& B, Matrix<T>& C)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    printf("Slave rank %d \n",rank);
    err = MPI_Recv(A.smat[rank], size2, MPI_INT, 0, ATAG, MPI_COMM_WORLD, &status);
    submatprintf(A.smat[rank], size, size, "Aslave", rank);
    CHECKMPI(err); 

    err = MPI_Recv(B.smat[rank], size2, MPI_INT, 0, BTAG, MPI_COMM_WORLD, &status);
    submatprintf(B.smat[rank], size, size, "Bslave", rank);
    CHECKMPI(err); 
    
    if(status.MPI_TAG == DIETAG) return;

    work(A.smat[rank], B.smat[rank], C.smat[rank]);

    err = MPI_Isend(C.smat[rank], size2, MPI_INT, 0, CTAG, MPI_COMM_WORLD, &request);
    CHECKMPI(err); 
}


