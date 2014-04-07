#include <array>
#include <stdio.h>
#include <iostream>
#include <string>

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
    }

    T& operator()(int i, int j) {return mat[j + i * cols];}

    Matrix<T>& operator=(const Matrix<T>& rhs)
    {
        if(mat != NULL)
            delete [] mat;
        mat = NULL;
        mat = new T[rhs.rows*rhs.cols];
        for(int i=0; i<rhs.rows; i++)
            for(int j=0; j< rhs.cols; j++)
                mat[j + i * cols] = rhs.mat[j + i * cols];
        return *this;
    }
    
    //T* submat(int row)
    //{
     //   return &this->mat[row*cols];
    //}

    void partition(int size, int N, int rank, bool transpose = false)
    {
        // size : size of the sub matrix
        // N : number of sub matrices per row or col
        // transpose : transpose the positions of submatrices

        int row, col, i, j, map=0;
        if (transpose)
            for(row=0; row<N; row+=size)
                for(col=row+size; col<N; col+=size)
                    for(i=0; i<size; ++i)
                        for (j=0; j<size; ++j)
                        {
                            int r1 = row+i; // original row index
                            int r2 = col+i; // transposed row index
                            int c1 = col+j; // original col index
                            int c2 = row+j; // transposed col index

                            // compute matrix indices
                            int idx1 = r1*N+c1;
                            int idx2 = r2*N+c2;

                            // swap indicies
                            T temp = mat[idx1];
                            mat[idx1] = mat[idx2];
                            mat[idx2] = temp;
                        }

        // reorganize mat
        cmat = new T[row*cols];
        for (row=0; row<N; row+=size)   // submatrix "row"
            for (col=0; col<N; col+=size) // submatrix "col"
                for (i=row; i<row+size; ++i) // full matrix row
                    for (j=col; j<col+size; ++j) // full matrix col
                        cmat[map++] = mat[i*N+j];
    }

    void clear()
    {
        int i;
        for(i=0;i<rows*cols;i++)
            mat[i] = 0;
             
    }
    
};

template <typename T>
void matmult(Matrix<T> A, Matrix<T>B, Matrix<T> C)
{
    int i,j,k;
    for(i=0;i<A.rows;i++)
        for(j=0;j<B.cols;j++)
            for(k=0;k<A.cols;k++) 
                C(i,j)+=A(i,k)*B(k,j);
}

template <typename T>
void matprintf(Matrix<T> A){
	int i,j,m,n;
    m = A.rows;
    n = A.cols;
	for(i=0;i<m;i++)
    {
        printf("\t\t");
        for(j=0;j<n;j++)    
            printf("%d ", A(i,j));
        printf("\n");
    }
}

template <typename T>
void submatprintf(T* A, int rows, int cols, std::string s = "", int rank = 0){
	int i,j;
	for(i=0;i<rows;i++)
    {
        printf("submat in %d %s\t ", rank, s.c_str());
        for(j=0;j<cols;j++)    
            printf("%d ", A[i * cols + j]);
        printf("\n");
    }
}
