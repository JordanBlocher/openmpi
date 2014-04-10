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
#define CTAG 3
#define DIETAG 4
#define MAXNUMPROCS 20

template <typename T>
class Matrix
{
 public:
    T* mat;
    T* cmat; // contiguous row-index matrix
    int rows;
    int cols;

    Matrix(int r, int c)
    {
       rows = r;
       cols = c;
       mat = new T[cols*rows]; 
       cmat = new T[cols*rows]; 
    }

    T& operator()(int i, int j) {return mat[j + i * cols];}
    
    Matrix<T>& operator=(const Matrix<T>& rhs)
    {
        if(mat != NULL)
            delete [] mat;
        if(cmat != NULL)
            delete [] cmat;
        mat = new T[rhs.rows*rhs.cols];
        cmat = new T[rhs.rows*rhs.cols];
        for(int i=0; i<rhs.rows*rhs.cols; ++i)
        {
            mat[i] = rhs.mat[i];
            cmat[i] = rhs.cmat[i];
        }
        return *this;
    }

    void shift(int size, int N, bool colshift = false, bool align = false)
    {
        // size : size of the sub matrix
        // N : number of sub matrices per row or col
        // transpose : transpose the positions of submatrices
        
        int row, col, i, r1, r2;
        T* temp = new T[rows*cols];
        int bmSize = N/size;
        int size2 = size*size;
        int bmSize2 = bmSize*bmSize;

        for(row=0; row<bmSize; ++row)
            for(col=0; col<size2; ++col)
                for(i=0; i<bmSize; ++i)
                {
                    if(colshift)
                    {
                        r1 = i*bmSize+row;
                        if(align)
                            r2 = ((i+row)*bmSize)%bmSize2+row;
                        else
                            r2 = ((i+1)*bmSize)%bmSize2+row;
                    }
                    else
                    {
                        r1 = row*bmSize+i;
                        if(align)
                            r2 = row*bmSize+(i+row)%bmSize;
                        else
                            r2 = row*bmSize+(i+1)%bmSize;
                    }
                    temp[r1*size2+col] = cmat[r2*size2+col];
                }
        
        delete [] cmat;
        cmat = temp;
   }

    void setRowIndexing(int size, int N, bool mapBack=false)
    {
        int row, col, i, j, map=0;

        // reorganize mat
        for (row=0; row<N; row+=size)   // submatrix "row"
            for (col=0; col<N; col+=size) // submatrix "col"
                for (i=row; i<row+size; ++i) // full matrix row
                    for (j=col; j<col+size; ++j) // full matrix col
                    {
                        if (mapBack)
                            cmat[map++] = mat[i*N+j];
                        else
                            mat[i*N+j] = cmat[map++];
                    }
    }

    void clear()
    {
        int i;
        for(i=0;i<rows*cols;i++)
        {
            mat[i] = 0;
            if(cmat != NULL)
                cmat[i] = 0;
        }
             
    }
    
};

template <typename T>
void pmatprintf(Matrix<T> A, int size, std::string s)
{
    A.setRowIndexing(size, A.rows, false);
	int i,j,m,n;
    m = A.rows < 15 ? A.rows : 15;
    n = A.cols < 15 ? A.cols : 15;
	for(i=0;i<m;i++)
    {
        if(i%size==0) printf("\t\t%*s\n", size*size, "-");
        printf("%s\t\t", s.c_str());
        for(j=0;j<n;j++)    
        {
            if(j%size==0) printf(" | ");
            printf("%d ", A(i,j));
        }
        printf("\n");
    }
}

template <typename T>
void work(T* A, T* B, T* C, int size)
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


int main (int argc, char *argv[]) 
{
    if(argc<2)
    {
        printf("Usage: ./bin/matrix <N> \n");
        return 0;
    }
    
    int N, numprocs, size, size2, bmSize;
    int rank, err=0;
    int i,j;
    N = atoi(argv[1]);
  
    std::chrono::high_resolution_clock::time_point startclock, stopclock;
    std::chrono::duration<double> time_span;
 
    // Init MPI & get numprocs and rank
    MPI_Init(&argc,&argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    CHECKMPI(err); 
    
    if(N % (int)sqrt(numprocs) != 0 || numprocs < 4 || fmod(sqrt(numprocs), 1) > 10e-6)
        return 0;

    size = N / (int)sqrt(numprocs);
    size2 = size*size;
    bmSize = N/size;
    Matrix<int> A(N,N);
    Matrix<int> B(N,N);
    Matrix<int> C(N,N);
    C.clear();

    if(rank == 0)
    {
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

        //pmatprintf(A, size, "A");
        //pmatprintf(B, size, "B");
        
        startclock = std::chrono::high_resolution_clock::now();

        A.shift(size, N, false, true);
        B.shift(size, N, true, true);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0; i<bmSize; ++i)
    {
        // Scatter B to processes
        err = MPI_Scatter(B.cmat, size2, MPI_INT, &(B.cmat[rank*size2]), size2, MPI_INT, 0, MPI_COMM_WORLD);
        CHECKMPI(err); 
        
        // Scatter A to processes
        err = MPI_Scatter(A.cmat, size2, MPI_INT, &(A.cmat[rank*size2]), size2, MPI_INT, 0, MPI_COMM_WORLD);
        CHECKMPI(err); 

        // Perform work on matrix
        work(&(A.cmat[rank*size2]), &(B.cmat[rank*size2]), &(C.cmat[rank*size2]), size);

        // Gather back result to C
        err= MPI_Gather(&(C.cmat[rank*size2]), size2, MPI_INT, C.cmat, size2, MPI_INT, 0, MPI_COMM_WORLD);
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
        
        if (i!=(bmSize-1) && rank == 0)
        {
            A.shift(size, N);
            B.shift(size, N, true);
        }
    }

    if(rank == 0)
    {
        stopclock = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stopclock - startclock);

    //    std::cout<<" = C\n";
      //  pmatprintf(C, size, "C");

        printf("\n %d %d %f\n",numprocs, N, (double)time_span.count());
    }

    MPI_Finalize();
    return 0;
}


