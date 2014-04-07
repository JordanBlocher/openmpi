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

MPI_Request request[20];

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
void work(T*, double);

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

    Construct_MPI_Datatypes<double>(rows, cols);
   
    if(rank == 0)
    {
        master<double>(numprocs, img);
    }
    slave<double>(numprocs);

    std::cout<<"rank " <<rank <<" waiting at barrier\n";
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        end = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        printf("\nProgram exec time %dx%d image: %f seconds \n\n",rows, cols, (double)time_span.count());

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
    Matrix<T> buf(numprocs+1, cols+1);
    //MPI_Request request;
    std::stringstream ss;
    MPI_Status status;
    MPI_Request request[numprocs+1];
  
    int i, j, k, err=0, complete, index;
    
    // Initialize processes
    for(i=0; i<numprocs; ++i)
    {
        buf[i][0] = i;
        std::cout<<"giving row "<<buf[i][0] << " to rank "<<i<<std::endl;
        err = MPI_Isend(buf[i], 1, MPI_Vector, i, WORKTAG, MPI_COMM_WORLD, &request[i]);
        CHECKMPI(err); 
    }
    
    k = 0;
    // Assign processes more rows
    int stop = std::ceil((rows+1)/2);
    while(k + numprocs < stop)
    {
            err = MPI_Testany(numprocs, request, &index, &complete, &status);
            CHECKMPI(err); 
            if(complete && index != MPI_UNDEFINED)
            {
                err = MPI_Recv(buf[index], 1, MPI_Vector, index, WORKTAG, MPI_COMM_WORLD, &status);
                std::cout<<"got row " << buf[index][0]<<" from " <<index + 1<< " \n";
                CHECKMPI(err); 
                
                    // Add data
                for(j=0; j<cols; j++)
                {
                    img[buf[index][0]][j] = (buf[index][j+1] < 0.001) ? 0 : buf[index][j+1];
             //       std::cout<<img[k][j]<<" "<< img[i][j] <<" " << img[i][j] << " ";
                }
                buf[index][0] = k + numprocs - 1;
                std::cout<<"giving row "<<buf[index][0] << " to index "<<index<<std::endl;
                err = MPI_Isend(buf[index], 1 ,MPI_Vector, index, WORKTAG, MPI_COMM_WORLD, &request[index]);
                CHECKMPI(err); 
                k++;
            }
    }

    //std::cout<<"No more rows\n";
    // No more rows, wait for processes to finish.
    for(i=0; i<numprocs; ++i)
    {
        err = MPI_Recv(buf[i], 1, MPI_Vector, MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD, &status);
            std::cout<<"got row " << buf[i][0]<<" from " <<status.MPI_SOURCE<< " \n";
        CHECKMPI(err); 
        for(j=0; j<cols; j++)
           img[buf[i][0]][j] = (buf[i][j+1] < 0.001) ? 0 : buf[i][j+1];
    }

    // Exit all slaves.
    for(i=0; i<numprocs; ++i)
    {
        err = MPI_Send(0, 0, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
        CHECKMPI(err); 
    }

    return;
}

template <typename T>
void slave(int numprocs)
{
    //Buffers
    Matrix<T> buf(numprocs+1, cols+1);
    double y;

    MPI_Status status;
    MPI_Request request;
    int rank, err=0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //std::cout<<"in slave rank " <<rank <<"\n";

    while(1)
    {
        err = MPI_Irecv(buf[rank], 1, MPI_Vector, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
        CHECKMPI(err); 
        if(status.MPI_TAG == DIETAG){return;}

        y = -2.0 + buf[rank-1][0]*(4.0/rows);
        work(buf[rank-1], y);
    // for(int j=1; j<cols; j++)
      //      std::cout<<buf[rank-1][j]<<" "<<std::trunc(buf[rank-1][j]/2)<<" "<<std::trunc(buf[rank-1][j]/2)<<" ";
     //std::cout<<"\n";

        err = MPI_Isend(buf[rank], 1, MPI_Vector, 0, WORKTAG, MPI_COMM_WORLD, &request);
        CHECKMPI(err); 
    }
}

    template <typename T>
void work(T* arr, double y)
{
    Complex z, z_next;
    z.Re = 0;; z.Im = 0;
    double x, hue, znorm, nu;
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
            arr[i+1] = hue;
        }
        else arr[i+1] = 0;
        //std::cout<<arr[i+1]<<" ";
        //if((z.Re*z.Re + z.Im*z.Im) > 4)
            //std::cout<<"z in set"<<z.Re<< " "<<z.Im <<std::endl<<std::endl;
    }
    //std::cout<<"\n";
}

    template <typename T>
void Construct_MPI_Datatypes(int rows, int cols)
{
    // Contiguous memory vector
    MPI_Type_contiguous(cols, MPI_DOUBLE, &MPI_Vector);
    MPI_Type_commit(&MPI_Vector);

    // Contiguous memory matrix
    MPI_Type_contiguous(rows, MPI_Vector, &MPI_Matrix);
    MPI_Type_commit(&MPI_Matrix);

    return;
}
