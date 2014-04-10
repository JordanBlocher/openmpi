#include <stdio.h>
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
#include <iostream>

template <typename T>
struct Matrix
{
    T** mat;

    Matrix(int rows, int cols)
    {
       mat = new T*[cols]; 
       for(int i=0; i<cols; i++)
       {
           mat[i] = new T[rows];
           for(int j=0; j< rows; j++)
               mat[i][j] = 0.0;
       }
        
    }

    T* operator[](const size_t i) {return mat[i];};
};

template <typename T>
void matmult(int n, Matrix<T> A, Matrix<T>B, Matrix<T> C);
template <typename T>
void matprintf(int n, Matrix<T> A);
template <typename T>
void vecprintf(int n, T *u);
int N;

int main (int argc, char *argv[]) {
    if(argc<2)
    {
        printf("Usage: ./executables/mandelbrot <N> \n");
        return 0;
    }

    int i, j;
    N = atoi(argv[1]);
    Matrix<float> A(N, N);
    Matrix<float> B(N, N);
    Matrix<float> C(N, N);
    std::string filename;
    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<float> time_span;
 
    start = std::chrono::high_resolution_clock::now();

    // Initialize matrix A random, B identity
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            A[i][j] = rand() % 20 + 1;
            if(i==j)
                B[i][j] = 1;
        }
    }

    matprintf(N, A);
    printf("X\n");
    matprintf(N, B);
    matmult(N, A, B, C);
    printf(" = \n");
    matprintf(N, C);


    end = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<float>>(end - start);
    printf("\nProgram exec time %dx%d image: %f seconds \n\n",N, N, (float)time_span.count());

    return 0;
}


template <typename T>
void matmult(int n, Matrix<T> A, Matrix<T>B, Matrix<T> C){
    int i,j,k;
    for(i=0;i<n;i++) for(j=0;j<n;j++)
        for(k=0;k<n;k++) 
            C[i][j]+=A[i][k]*B[k][j];
}

template <typename T>
void matprintf(int n, Matrix<T> A){
	int i,m;
	if(n>5) m=5;
	else m=n;
	for(i=0;i<m;i++) vecprintf<T>(n,A[i]);
    if(n!=m) printf("\n .\n. \n.");
}

template <typename T>
void vecprintf(int n, T *u)
{
    int i,m;
	if(n<5) m=n;
	else m=5;
    for(i=0;i<m;i++) printf("%5.5f ",u[i]);
	if(n!=m) printf("...");
    printf("\n");
}
