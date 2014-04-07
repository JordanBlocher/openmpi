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

#define WORKTAG 1
#define DIETAG 2
#define MAXDATASIZE 500
#define MAXNUMPROCS 20
#define ITER 256

typedef std::array<int,3> color;
std::array<color, 5> Palette = {{ color {{1, 20, 30}}, color {{1, 20, 10}}, color {{20, 0, 30}}, color {{10, 0, 10}}}};


//MPI Datatypes
MPI_Datatype MPI_Vector;
MPI_Datatype MPI_Matrix;
MPI_Datatype MPI_Complex;


template <typename T>
struct Matrix
{
    T** mat;

    Matrix(int rows, int cols)
    {
       mat = new T*[rows]; 
       for(int i=0; i<rows; i++)
       {
           mat[i] = new T[cols];
           for(int j=0; j< cols; j++)
               mat[i][j] = 0.0;
       }
        
    }

    T* operator[](const size_t i) {return mat[i];};
};

struct Complex
{
    double Re;
    double Im;

    void operator = (const Complex& z) {Re = z.Re; Im = z.Im;};
};

struct Pixel
{
    int x;
    int y;
    double hue[3];
};

template <typename T>
void master(int, Matrix<T>&);
template <typename T>
void slave(int);
template <typename T>
void Construct_MPI_Datatypes(int, int);
template <typename T>
void work(T*, double, int);

int rows, cols;

int main (int argc, char *argv[]) {
    if(argc<3)
    {
      //  printf("Usage: ./executables/mandelbrot <rows> <cols> (filename)\n");
        return 0;
    }

    int rank, numprocs, err=0, i, j;
    std::vector<std::string> v;
    rows = atoi(argv[1]);
    cols = atoi(argv[2]);
    Matrix<double> img(rows, cols);
    std::string filename;
    if(argc > 3)
    {
        filename = argv[3];
        filename += ".ppm";
    }

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> time_span;
 
    start = std::chrono::high_resolution_clock::now();
    // Init MPI & get numprocs and rank
    MPI_Init(&argc,&argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if(err){
        fprintf(stderr,"Catastrophic MPI problem.\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    //int gridSize = std::min((int)std::ceil(rows/(numprocs-1)), 1000);
    //if(gridSize % rows != 0) gridSize = (rows % 1000 == 0) ? 1000 : 500;
    int gridSize = 500;
    Construct_MPI_Datatypes<double>(gridSize, cols);
   
    if(rank == 0)
    {
        master<double>(numprocs, img);
    }
    if(rank != 0){
        slave<double>(numprocs);
    }

    std::cout<<"rank " <<rank <<" waiting at barrier\n";
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        end = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        printf("\nProgram exec time %dx%d image %d processors: %f seconds \n\n",rows, cols, numprocs, (double)time_span.count());

        if(argc > 3)
        {
            std::ofstream fout;
            fout.open(filename, std::ofstream::out | std::ofstream::binary);
            fout << "P3\n # Mandelbrot\n " << rows << " " <<cols << "\n 255\n";
             double val;
            int r, g, b;
            for(i=0; i<=std::floor(rows/2.0); i++)
            {
                for(j=0; j<cols; j++)
                {
                    val = std::abs(img[i][j]);
                    r = val * Palette[(int)val % 5][0]; 
                    g = val * Palette[(int)val % 5][1]; 
                    b = val * Palette[(int)val % 5][2]; 
                    fout << r << " " << g <<" " << b << " ";
                }
                fout <<"\n";
            }
            for(i=std::floor(rows/2.0); i>=0; i--)
            {
                for(j=0; j<cols; j++)
                {
                    val = std::abs(img[i][j]);
                    r = val * Palette[(int)val % 5][0]; 
                    g = val * Palette[(int)val % 5][1]; 
                    b = val * Palette[(int)val % 5][2]; 
                    fout << r << " " << g <<" " << b << " ";
                }
                fout <<"\n";
            }
            fout.close();
        }
    }

    // Free memory
    MPI_Type_free(&MPI_Vector);

    MPI_Finalize();
    return 0;
}

template <typename T>
void master(int numprocs, Matrix<T>& img)
{

    // Buffers
    //int gridSize = std::min((int)std::ceil(rows/(numprocs-1)), 1000);
    //std::cout<<"gridsize "<<gridSize<<"\n";
    //if(gridSize % rows != 0) gridSize = (rows % 1000 == 0) ? 1000 : 500;
    int gridSize = 500;
    Matrix<T> buf(numprocs, gridSize*cols+1);
    MPI_Request request[numprocs];
    MPI_Status status;
    std::stringstream ss;
    std::cout<<"gridsize "<<gridSize<<"\n";
  
    int stop = std::ceil((rows+1)/2.0);
    int i, j, k, err=0, cnt=0, idx, complete=1, rank;
    
    // Initialize processes
    for(i=0; i<numprocs - 1 && cnt < stop; ++i)
    {
        buf[i][0] = cnt;
        std::cout<<"giving row "<<buf[i][0] << " to " << cnt + gridSize <<" to rank "<<i+1<<std::endl;
        err = MPI_Isend(buf[i], 1, MPI_Vector, i+1, WORKTAG, MPI_COMM_WORLD, &request[i]);
        CHECKMPI(err); 
        cnt+=gridSize;
        idx = i+1;
    }
    if(idx < numprocs - 1)
    {
        // Exit useless slaves.
        for(i=idx; i<numprocs - 1; ++i)
        {
            std::cout<<"deleting "<<i+1<<"\n";
            err = MPI_Isend(0, 0, MPI_INT, i+1, DIETAG, MPI_COMM_WORLD, &request[i]);
            CHECKMPI(err); 
        }
    }
    numprocs = idx;
    std::cout<<"idx "<<idx<<"\n";
    // Assign processes more rows
    while(cnt <= stop)
    {
        err = MPI_Testany(numprocs, request, &rank, &complete, &status);
        CHECKMPI(err); 
        std::cout<<"(loop )complete= "<<complete<<" in slave "<<status.MPI_SOURCE<<"\n";
        if(complete && idx != MPI_UNDEFINED)
        {
            if(cnt >= stop) break;
            err = MPI_Recv(buf[rank], 1, MPI_Vector, rank, WORKTAG, MPI_COMM_WORLD, &status);
            std::cout<<"got row " << buf[i][0]<< " to " << buf[i][0]+gridSize <<" from rank " <<status.MPI_SOURCE<< " \n";
            CHECKMPI(err); 
            // Add data
            idx = 0;
            for(k=0; k<gridSize; k++)
            for(j=0; j<cols; j++)
            {
                img[buf[i][0]+k][j] = (buf[i][idx+1] < 0.001) ? 0 : buf[i][idx+1];
                idx++;
       //         std::cout<<buf[i-1][j+1]<<" "<< buf[i-1][j+1] <<" " << buf[i-1][j+1] << " ";
            }
            std::cout<<"\n";
            buf[i][0] = cnt;
            std::cout<<"giving row "<<buf[i][0] << " to " << buf[i][0] + gridSize << " to rank "<<status.MPI_SOURCE<<std::endl;
            err = MPI_Isend(buf[rank], 1 ,MPI_Vector, rank, WORKTAG, MPI_COMM_WORLD, &request[rank]);
            CHECKMPI(err); 
            cnt+=gridSize;
        }
        complete = 1;
    }

    std::cout<<"No more rows\n";
    complete = 1;
    // No more rows, wait for processes to finish.
    for(i=0; i<numprocs; ++i)
    {
        err = MPI_Testany(numprocs, request, &rank, &complete, &status);
        CHECKMPI(err); 
        std::cout<<"complete= "<<complete<<" in slave "<<status.MPI_SOURCE<<"\n";
        if(complete && idx != MPI_UNDEFINED)
        {
            err = MPI_Recv(buf[i], 1, MPI_Vector, i+1, WORKTAG, MPI_COMM_WORLD, &status);
                std::cout<<"got row " << buf[i][0]<<" to "<<buf[i][0] + gridSize<<" from " <<status.MPI_SOURCE<< " \n";
            CHECKMPI(err); 
            idx = 0;
            for(k=buf[i][0]; k<buf[i][0]+gridSize; k++)
            {
            for(j=0; j<cols; j++)
            {
               img[k][j] = (buf[i][idx+1] < 0.001) ? 0 : buf[i][idx+1];
               idx++;
            }
            }
        }
        complete = 1;
    }

    // Exit all slaves.
    for(i=0; i<numprocs; ++i)
    {
        std::cout<<"Exiting slaves \n";
        err = MPI_Isend(0, 0, MPI_INT, i+1, DIETAG, MPI_COMM_WORLD, &request[i]);
        CHECKMPI(err); 
    }

    return;
}

template <typename T>
void slave(int numprocs)
{
    //Buffers
    //int gridSize = std::min((int)std::ceil(rows/(numprocs-1)), 1000);
    //if(gridSize % rows != 0) gridSize = (rows % 1000 == 0) ? 1000 : 500;
    int gridSize = 500;
    Matrix<T> buf(numprocs, gridSize*cols+1);

    MPI_Status status;
    MPI_Request request;
    int rank, err=0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout<<"in slave rank " <<rank <<"\n";

    while(1)
    {
        err = MPI_Recv(buf[rank-1], 1, MPI_Vector, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        CHECKMPI(err); 
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(status.MPI_TAG == DIETAG){return;}

        work(buf[rank-1], gridSize, buf[rank-1][0]);
     //for(int j=1; j<cols; j++)
       //     std::cout<<buf[rank-1][j]<<" "<<std::trunc(buf[rank-1][j]/2)<<" "<<std::trunc(buf[rank-1][j]/2)<<" ";
     //std::cout<<"\n";

        err = MPI_Isend(buf[rank-1], 1, MPI_Vector, 0, WORKTAG, MPI_COMM_WORLD, &request);
        CHECKMPI(err); 
    }
}

    template <typename T>
void work(T* arr, double gridSize, int index)
{
    Complex z, z_next;
    z.Re = 0;; z.Im = 0;
    double x, y, hue, znorm, nu;
    int cnt = 0;
    for(int j=index; j<index+gridSize; j++)
    {
        y = -2.0 + j*(4.0/rows);
        for(int i=0; i<cols; i++)
        {
            x = -1.5 + i*(4.0/cols); 
            z.Re = 0;; z.Im = 0;
            for(int j=0; j<ITER && (z.Re*z.Re + z.Im*z.Im) <= 4 ; j++)
            {
                /* z^2 = (a+bi)(a+bi) = a^2 + 2abi - b^2 */
                z_next.Re = (z.Re*z.Re)-(z.Im*z.Im) + x;
                z_next.Im = 2*z.Re*z.Im + y;
                hue = j;
                z = z_next;
            }
            znorm = z.Re*z.Re + z.Im*z.Im;
            if(znorm >= 4)
            {
                nu = log( log(znorm) / log(2) ) / log(2); 
                hue = hue + 1 - nu;
                arr[cnt+1] = hue;
            }
            else arr[cnt+1] = 0;
            cnt++;
        }
    }
}

    template <typename T>
void Construct_MPI_Datatypes(int rows, int cols)
{
    // Contiguous memory vector
    MPI_Type_contiguous(cols*rows+1, MPI_DOUBLE, &MPI_Vector);
    MPI_Type_commit(&MPI_Vector);

    return;
}
