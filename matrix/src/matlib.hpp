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

