#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "string.h"
#include <time.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <sstream>

#define WORKTAG 1
#define DIETAG 2
#define MAXDATASIZE 20

//MPI Datatypes
MPI_Datatype MPI_String;
MPI_Datatype MPI_Vector;

struct String
{
    char c[MAXDATASIZE];
};

void master(int);
void slave(int);
void Construct_MPI_String(int);
template <typename T>
void Construct_MPI_Vector(T, int);

int main (int argc, char *argv[]) {
    if(argc<2)
    {
        printf("Usage: ./executables/pingpong <numprocs> \n");
        return 0;
    }

    int rank, numprocs, err=0;
    numprocs = atoi(argv[1]);
    clock_t start, end;
    std::vector<std::string> v;

    start = MPI_Wtime();

    // Init MPI & get numprocs and rank
    MPI_Init(&argc,&argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if(err){
        fprintf(stderr,"Catastrophic MPI problem.\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    Construct_MPI_String(MAXDATASIZE);
    //Construct_MPI_Vector(MPI_String, 1);

    if(rank == 0)
    {
        master(numprocs);
    }
    if(rank != 0){
        slave(numprocs);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Free memory
    MPI_Type_free(&MPI_String);
    //MPI_Type_free(&MPI_Vector);

    MPI_Finalize();

    end = MPI_Wtime();
    if(rank == 0){
        printf("\nAverage message time (one-way): %f seconds \n\n",(double)(end - start)/(double)(2.0*numprocs));}

    return 0;
}

void master(int numprocs)
{

    // Buffers
    std::vector<std::string> sendbuf;
    std::vector<String> recvbuf;
    //MPI_Request request;
    MPI_Status status;
    std::stringstream ss;
    String s;

    int i, err=0;
    std::vector<clock_t> start, end;

    for(i=0; i<numprocs-1; ++i)
    {
        ss << i;
        sendbuf.push_back(ss.str() + ": Ping -> ");
        ss.str("");
        recvbuf.push_back(s);
        start.push_back(MPI_Wtime());

        err = MPI_Send(&sendbuf[i][0], 1, MPI_String, i+1, WORKTAG, MPI_COMM_WORLD);
        if(err){
            fprintf(stderr,"Failed to send.\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }

    // Wait for processes to finish.
    for(i=0; i<numprocs-1; ++i)
    {
        err = MPI_Recv(&recvbuf[i], 1, MPI_String, MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD, &status);
        printf("%s\n", recvbuf[i].c);
        if(err){
            fprintf(stderr,"Failed to recieve.\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }

    // Exit all slaves.
    for(i=0; i<numprocs-1; ++i)
    {
        err = MPI_Send(0, 0, MPI_INT, i+1, DIETAG, MPI_COMM_WORLD);
        if(err){
            fprintf(stderr,"Catastrophic MPI problem.\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }

    return;
}

void slave(int numprocs)
{
    //Buffers
    std::vector<std::string> sendbuf;
    std::vector<String> recvbuf;
    String s;
    std::vector<clock_t> start, end;

    
    MPI_Status status;
    //MPI_Request request;
    int rank, err=0;

    for(int i=0; i< numprocs-1;i++)
    {
        sendbuf.push_back("");
        recvbuf.push_back(s);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    while(1)
    {
        err = MPI_Recv(&recvbuf[rank-1], 1, MPI_String, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if(err){
            fprintf(stderr,"Failed to receive.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(status.MPI_TAG == DIETAG){return;}

        sendbuf[rank-1] = std::string(recvbuf[rank-1].c) + "Pong\0";
        err = MPI_Send(&sendbuf[rank-1][0], 1, MPI_String, 0, WORKTAG, MPI_COMM_WORLD);
        if(err){
            fprintf(stderr,"Failed to send.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
}

void Construct_MPI_String(int count)
{
    // Contiguous memory string
    MPI_Type_vector(count, MAXDATASIZE, 0, MPI_CHAR, &MPI_String);
    MPI_Type_commit(&MPI_String);

    return;
}

    template <typename T>
void Construct_MPI_Vector(T type, int count)
{
    // Contiguous memory vector
    MPI_Type_indexed(count, type, &MPI_Vector);
    MPI_Type_commit(&MPI_Vector);

    return;
}
