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
                        // if !align
                        // row = 0 0 0 1 1 1 2 2 2
                        // i   = 0 1 2 0 1 2 0 1 2
                        // r1  = 0 3 6 1 4 7 2 5 8
                        // r2  = 3 6 0 4 7 1 5 8 2
                        // else
                        // row = 0 0 0 1 1 1 2 2 2
                        // i   = 0 1 2 0 1 2 0 1 2
                        // r1  = 0 3 6 1 4 7 2 5 8
                        // r2  = 0 3 6 4 7 1 8 2 5
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
        
        /*
        for(row=0; row<N; row+=size)
            for(col=0; col<N; col+=size)
                for(i=0; i<size; ++i)
                    for (j=0; j<size; ++j)
                    {
                        int r1 = row+i; // original row index
                        int r2 = row+i; // shifted row index
                        int c1 = col+j; // original col index
                        int c2 = col+j; // shifted col index
                        if(colshift)
                        {
                            cshift = align ? col : size;
                            r2 = (r2 - cshift) < 0 ? N + (r2 - cshift) : r2 - cshift;
                        }
                        else
                        {
                            rshift = align ? row : size;
                            c2 = (c2 - rshift) < 0 ? N + (c2 - rshift) : c2 - rshift;
                        }

                        // compute matrix indices
                        int idx1 = r1*N+c1;
                        int idx2 = r2*N+c2;

                        temp[idx2] = mat[idx1];
                    }
        */
 
        //setRowIndexing(size, N);
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
void pmatprintf(Matrix<T> A, int size, std::string s = "")
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

