#include <mpi.h>
#include <cstdio>
#include <cstdlib>

using std::cout;
using std::endl;


#define MASTER 0


typedef struct
{
  int my_rank;
  int total_nodes;
  int 
} JobInfo;



int master( JobInfo& info, string& input_name );


int main( int argc, char** argv )
{
    // variables
    JobInfo info;

    // do MPI initialization
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &(info.total_nodes) );
  MPI_Comm_rank( MPI_COMM_WORLD, &(info.my_rank) );

    

    // case: this processor is the master
    if( info.my_rank == MASTER )
    {
        // run the master's code
        master(         );
    }
    // case: this node is not the master
    else
    {

    }

    // finalize the program
    MPI_Finalize( MPI_COMM_WORLD );

    // return the program status indicator
    return 0;
}


int master( JobInfo& info, string& input_name )
{
    // read in the file data

    // 
}


