#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstring>
#include <sstream>
#include <mpi.h>
#include <chrono>
#include "my_stopwatch.h"

using std::cout;
using std::endl;

// standard constants
#define STD_STR_LEN 40

// node designations
#define MASTER 0

// timers
#define TOTAL_CLOCK 0

// n-body specific constants
#define TIME_DIFFERENCE 0.0003
#define G 6.67384e-11
#define DT 100 // Big time step so to see movement

#define ARBITRARY_TAG 1


// global variables
MPI_Datatype MPI_Body;

typedef struct
{
    int my_rank;
    int total_nodes;
    char input_file_name[STD_STR_LEN];
    int num_simulation_steps;
    float time_difference;
    int num_bodies;
    int* start_indices;
    int* shares_of_bodies;
} JobInfo;

struct Body
{
    int mass;
    float position[3];
    float velocity[3];    
};

struct Force
{
    float magnitude[3];  // each magnitude corresponds to a cartesian direction
};

struct Node;

union Cell
{
    Body *body;
    Node *nodes[8];
};

struct Node
{
    int type;
    float mass;
    float position[3]; // center of mass
    Cell cell;
};


void CreateMPIDatatype();

int processCommandLineArgs( int argc, char** argv, JobInfo* job_info );

void computeJobShares( JobInfo* info );

int master( JobInfo* info );

int slave( JobInfo* info );

int readNbodyData( JobInfo* info, Body** universe );

void transmitProblemParameters( JobInfo* info );

void receiveProblemParameters( JobInfo* info );

void allocateArrays( JobInfo* info, Body** universe, Body** my_bodies,
                     Body** results );

void buildOctreeFromData(  /* PUT STUFF HERE!!!*/  );

void loadMyBodies( Body* my_bodies, Body* universe, JobInfo* info );

void bodyCopy( Body* destination, Body* source );

void performComputations( Body* results, Body* my_bodies, Body* universe,
                          JobInfo* info );

void addTo3dForce( Force* resulting_force, Body* reference_body,
                   Body* other_body );

void createResultingBody( Body* result_body, Body* original_body,
                          Force* cumulative_force );

void findNewVelocity( Body* result_body, Body* original_body,
                      Force* cumulative_force );

void findNewPosition( Body* result_body, Body* original_body );

void outputPerformanceSummary( JobInfo* info, float total_simulation_time );

void printUniverse( Body* universe );

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
    // variables
    JobInfo info;

    // do MPI initialization
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &(info.total_nodes) );
    MPI_Comm_rank( MPI_COMM_WORLD, &(info.my_rank) );
    CreateMPIDatatype();

    // process the command line arguments
    processCommandLineArgs( argc, argv, &info );

// this will be unnecessary once octree is in place
computeJobShares( &info );





    // case: this processor is the master
    if( info.my_rank == MASTER )
    {
        // run the master's code
        master( &info );
    }
    // case: this node is not the master
    else
    {
        // run the slave's code
        slave( &info );
    }

printf( "MADE IT TO LINE %d\n", __LINE__ );
fflush( stdout );

    // finalize the program
    MPI_Finalize();

    // return the program status indicator
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CreateMPIDatatype()
{

/*
    int blocks = 3;
    
    MPI_Type_struct( 3,
// count  // array of blocklengths // array of displacements  // array of types  // newtype )
  */

    // Array of datatypes
    MPI_Datatype type[3] = {  MPI_INT, MPI_FLOAT, MPI_FLOAT };

    // Number of each type
    int block_len[3] = { 1, 3, 3};
    
    // Displacements
    MPI_Aint disp[3] = { 0, sizeof(MPI_INT),
                         (sizeof(MPI_INT) + (3 * sizeof(MPI_FLOAT))) };

    // MPI Body struct
    MPI_Type_create_struct(3, block_len, disp, type, &MPI_Body);

    MPI_Type_commit(&MPI_Body);
}

int processCommandLineArgs( int argc, char** argv, JobInfo* job_info )
{
    // variables
    int arg_success = false;

    // case: there is a sufficient number of command line arguments
    if( argc >= 2 )
    {
        // incorporate the arguments into the program
        strcpy( job_info->input_file_name, argv[1] );

        // indicate that arguments were successfully obtained
        arg_success = true;
    }

    // return the status of reading in of the command line arguments
    return arg_success;
}


void computeJobShares( JobInfo* info )
{
  // variables
  info->start_indices = NULL;
  info->shares_of_bodies = NULL;
  int minimum_share = 0;
  int bigger_share = 0;
  int remainder = 0;
  int i = 0;

  // allocate the row information arrays
  info->start_indices = (int*) calloc( info->total_nodes, sizeof( int ) );
  info->shares_of_bodies = (int*) calloc( info->total_nodes, sizeof( int ) );

  // compute the minimum share and the bigger share of the rows
  minimum_share = info->num_bodies / info->total_nodes;
  remainder = info->num_bodies % info->total_nodes;
  bigger_share = minimum_share + 1;

  // distribute the minimum share and the remainder among the low rank nodes
  for( i = 0; i < remainder; i++ )
  {
    // allocate the start row and row shares
    info->start_indices[i] = i * bigger_share;
    info->shares_of_bodies[i] = bigger_share;
  }

  // finish distributing the starting rows and shares
  for( /* i = whatever it was from last loop */; i < info->total_nodes; i++ )
  {
    // allocate the start row and row shares
    info->start_indices[i] = ( i * minimum_share ) + remainder;
    info->shares_of_bodies[i] = minimum_share;
  }

  // return the job info by reference
}


int master( JobInfo* info )
{
    // variables
    int job_success = true;
    int i = 0;
    Body* universe = NULL;
    Body* my_bodies = NULL;
    Body* results = NULL;
    float total_simulation_time = 0;

printf( "MADE IT TO LINE %d\n", __LINE__ );
fflush( stdout );

    // read in the file data
    readNbodyData( info, &universe );

printf( "MADE IT TO LINE %d\n", __LINE__ );
fflush( stdout );


    // allocate the arrays appropriately
    my_bodies = ( Body* ) calloc( info->shares_of_bodies[ info->my_rank ],
                                  sizeof( Body ) );
    results = ( Body* ) calloc( info->shares_of_bodies[ info->my_rank ],
                                sizeof( Body ) );

printf( "MADE IT TO LINE %d\n", __LINE__ );
fflush( stdout );

    // transmit the pertinent problem data to all processes
    transmitProblemParameters( info );

printf( "MADE IT TO LINE %d\n", __LINE__ );
fflush( stdout );

printf( "num bodies %d\n", info->num_bodies );
fflush( stdout );

    // begin the timing operation
    MPI_Barrier( MPI_COMM_WORLD );

    // scatter/broadcast the initial data
    MPI_Bcast( universe,           // send buffer
               info->num_bodies,   // number of items
               MPI_Body,           // send type
               MASTER,             // root
               MPI_COMM_WORLD );   // communicator

printf( "MADE IT TO LINE %d\n", __LINE__ );
fflush( stdout );


    // do n-body computations until the simulation has run for enough time
    for( i = 0; i < info->num_simulation_steps; i++ )
    {
        // re-build the octree
        buildOctreeFromData(  /* STUFF WILL GO HERE */  );

        // load this node's set of bodies for the computations
        loadMyBodies( my_bodies, universe, info );

        // do the computations
        performComputations( results, my_bodies, universe, info );

        // collect and share results from all processors
        MPI_Allgatherv(
            results,                            // send buffer
            info->shares_of_bodies[ info->my_rank ], // send count
            MPI_Body,                           // send type
            universe,                           // receive buffer
            info->shares_of_bodies,                  // receive counts
            info->start_indices,                // displacements
            MPI_Body,                           // receive type
            MPI_COMM_WORLD );                   // communicator
    }


printf( "MADE IT TO LINE %d\n", __LINE__ );
fflush( stdout );

    // stop the timing
    MPI_Barrier( MPI_COMM_WORLD );
    total_simulation_time = stopwatch( 'x', TOTAL_CLOCK );

printf( "MADE IT TO LINE %d\n", __LINE__ );
fflush( stdout );

    // output the performance summary
    outputPerformanceSummary( info, total_simulation_time );

    // return the success status of the job
    return job_success;
}


int slave( JobInfo* info )
{
    // variables
    int job_success = true;
    int i = 0;
    Body* universe = NULL;
    Body* my_bodies = NULL;
    Body* results = NULL;

printf( "SLAVE LINE %d\n", __LINE__ );
fflush( stdout );

    // collect the problem parameters
    receiveProblemParameters( info );

printf( "SLAVE LINE %d\n", __LINE__ );
fflush( stdout );

    // prep the arrays for the problem
    allocateArrays( info, &universe, &my_bodies, &results );

printf( "SLAVE LINE %d\n", __LINE__ );
fflush( stdout );

    // synchronize with other processes to properly start timing
    MPI_Barrier( MPI_COMM_WORLD );

printf( "SLAVE LINE %d\n", __LINE__ );
fflush( stdout );

printf( "num bodies %d\n", info->num_bodies );
fflush( stdout );


    // receive a data set from the master
    MPI_Bcast( universe,           // send buffer
               info->num_bodies,   // number of items
               MPI_Body,           // send type
               MASTER,             // root
               MPI_COMM_WORLD );   // communicator

printf( "SLAVE LINE %d\n", __LINE__ );
fflush( stdout );

    // do the n-body computations for the length of the simulation
    for( i = 0; i < info->num_simulation_steps; i++ )
    {
        // load this node's set of bodies for the computations
        loadMyBodies( my_bodies, universe, info );

        // do the computations
        performComputations( results, my_bodies, universe, info );

        // collect and share results from all processors
        MPI_Allgatherv(
            results,                            // send buffer
            info->shares_of_bodies[ info->my_rank ], // send count
            MPI_Body,                           // send type
            universe,                           // receive buffer
            info->shares_of_bodies,                  // receive counts
            info->start_indices,                // displacements
            MPI_Body,                           // receive type
            MPI_COMM_WORLD );                   // communicator
    }

    // synchronize with other processes to properly end timing
    MPI_Barrier( MPI_COMM_WORLD );

    // return the job success status
    return job_success;
}


int readNbodyData( JobInfo* info, Body** universe )
{
    // variables 
    int read_success = false;
    FILE* file = NULL;
    int i = 0, j = 0;
    Body body; // temp for file reading

    // open the file
    file = fopen( info->input_file_name , "r");

    // case: the file opened successfully
    if( file != NULL )
    {
        // read in the problem parameters
        fscanf( file, "%d %d %f", &(info->num_bodies),
                &(info->num_simulation_steps), &(info->time_difference) );

        // allocate the array to contain pointers to all of the bodies
        *universe = ( Body* ) calloc( info->num_bodies, sizeof( Body* ) );

        // read in all body information
        for( i = 0; i < info->num_bodies; ++i )
        {
            // read in the body's position as a spacial coordinate triple
            fscanf( file, "%f %f %f", &(body.position[0]),
                    &(body.position[1]), &(body.position[2]) );

            // read in the initial velocity for the body
            fscanf( file, "%f %f %f", &(body.velocity[0]),
                    &(body.velocity[1]), &(body.velocity[2]) );

            // read in the mass value
            fscanf( file, "%d", &(body.mass) );

            // store the newly read in body
            bodyCopy( &(body), &((*universe)[i]) ); 
        }

        // indicate that the file was successfully read
        read_success = true;
    }
    // case: the file did not open
    else
    {
        // give a nasty message
        printf( "BAD FILE NAME GUY!\r\n" );
    }

    // return file read success, return the universe array by reference
    return read_success;
}


void transmitProblemParameters( JobInfo* info )
{
    // variables
    int i = 0;

    // send each of the pertinent items to each processor
    for( i = 1; i < info->total_nodes; i++ )
    {
        // send the items
        MPI_Send( &(info->num_bodies), 1, MPI_INT, i, ARBITRARY_TAG,
                  MPI_COMM_WORLD );
        MPI_Send( &(info->num_simulation_steps), 1, MPI_INT, i, ARBITRARY_TAG,
                  MPI_COMM_WORLD );
        MPI_Send( &(info->time_difference), 1, MPI_FLOAT, i, ARBITRARY_TAG,
                  MPI_COMM_WORLD );
    }
}


void receiveProblemParameters( JobInfo* info )
{
    // variables
    MPI_Status status;

    // receive the items
    MPI_Recv( &(info->num_bodies), 1, MPI_INT, MASTER, MPI_ANY_TAG,
              MPI_COMM_WORLD, &status );
    MPI_Recv( &(info->num_simulation_steps), 1, MPI_INT, MASTER, MPI_ANY_TAG,
              MPI_COMM_WORLD, &status );
    MPI_Recv( &(info->time_difference), 1, MPI_FLOAT, MASTER, MPI_ANY_TAG,
              MPI_COMM_WORLD, &status );    
}


void allocateArrays( JobInfo* info, Body** universe, Body** my_bodies,
                     Body** results )
{
    // allocate the arrays as necessary
    *universe = ( Body* ) calloc( info->num_bodies, sizeof( Body ) );
    *my_bodies = ( Body* ) calloc( info->shares_of_bodies[ info->my_rank ],
                                   sizeof( Body ) );
    *results = ( Body* ) calloc( info->shares_of_bodies[ info->my_rank ],
                                 sizeof( Body ) );
}


void buildOctreeFromData(  /* STUFF WILL GO HERE */  )
{


        /* STUFF WILL GO HERE */


}


void loadMyBodies( Body* my_bodies, Body* universe, JobInfo* info )
{
    // variables
    int i = 0;

    // simply copy the bodies from the universe into this node's body set
    // using appropriate offsets
    for( i = 0; i < info->shares_of_bodies[ info->my_rank ]; i++ )
    {
        bodyCopy( &(my_bodies[i]),
                  &(universe[ i + info->start_indices[ info->my_rank ] ]) );
    }

    // The array of bodies for this node is returned by reference
}


void bodyCopy( Body* destination, Body* source )
{
    // variables
    int i = 0;

    // simply assign values from source to all members of destination 
    destination->mass = source->mass;
    for( i = 0; i < 3; i++ )
    {
        destination->position[i] = source->position[i];
    }
 
    for( i = 0; i < 3; i++ )
    {
        destination->velocity[i] = source->velocity[i]; 
    }
}


void performComputations( Body* results, Body* my_bodies, Body* universe,
                          JobInfo* info )
{
    // variables
    Force cumulative_force;
    int i = 0;
    int j = 0;

    // for each of the bodies
    for( i = 0; i < info->shares_of_bodies[ info->my_rank ]; i++ )
    {
        // case: the "reference" and "other" bodies are the same
        if( i == ( j - info->start_indices[ info->my_rank ] ) )
        {
            // skip this computation
            continue;
        }

        // reset the force vector
        for( j = 0; j < 3; j++ )
        {
            cumulative_force.magnitude[i] = 0;
        }

        // for every other body in the universe
        for( j = 0; j < info->num_bodies; j++ )
        {
            // account for the other body's force on this one
            addTo3dForce( &cumulative_force, &(my_bodies[i]), &(universe[j]) );
        }

        // we didn't check to see for the current body in the force computation
        // for performance reasons, subtract that factor back out
        // (this is a hacky way to accomplish the subtraction)
// DO THIS !!! There is a hack to get by for now ( the continue )

        // create and store the body state that results after the small time
        // difference due to the physics
        createResultingBody( &(results[i]), &(my_bodies[i]),
                             &cumulative_force );
    }

    // the results are returned by reference
}


void addTo3dForce( Force* resulting_force, Body* reference_body,
                   Body* other_body )
{
    // variables
    int i = 0;
    float radius_component = 0;
    float radius_squared = 0;

    // Note: r = sqrt( \sum{ ( d_2 - d_1 )^2 } ), but since r^2 is needed, we
    //       will not do the sqrt part 
    // find the distance between the two objects
    for( i = 0; i < 3; i++ )
    {
        // find a component of the radius (squared)
        radius_component = reference_body->position[i] -
                           other_body->position[i];
        radius_component = radius_component * radius_component;

        // add it to the r^2 part
        radius_squared += radius_component;
    }

    // for each dimension compute the force component
    for( i = 0; i < 3; i++ )
    {
        // apply F_x = G * { m_1 * m_2 } / { r^2 }
        resulting_force->magnitude[i] +=
            G * ( reference_body->mass * other_body->mass ) / radius_squared;
    }

    // the modified force exerted on the reference body is returned by reference 
}


void createResultingBody( Body* result_body, Body* original_body,
                          Force* cumulative_force )
{
    // variables
    int i = 0;

    // set up the masses of the result body
    result_body->mass = original_body->mass;

    // compute the new velocity
    findNewVelocity( result_body, original_body, cumulative_force );

    // compute the new position
    findNewPosition( result_body, original_body );

    // results are returned by reference
}


void findNewVelocity( Body* result_body, Body* original_body,
                      Force* cumulative_force )
{
    // variables
    int i = 0;

    // adjust the velocity components by applying the accelerations resulting
    // from the gravitational forces
    for( i = 0; i < 3; i++ )
    {
        // modify the velocity by the acceleration over the small time
        // difference v_f = v_i + ( F / m ) * dt
        result_body->velocity[i] =
            ( original_body->velocity[i] +
              ( ( cumulative_force->magnitude[i] / original_body->mass ) *
                TIME_DIFFERENCE ) );
    }

    // the result body is returned by reference
}


void findNewPosition( Body* result_body, Body* original_body )
{
    // variables
    int i = 0;

    // adjust the position components by adding the accelerations resulting
    // from the gravitational forces
    for( i = 0; i < 3; i++ )
    {
        // modify the position by the velocity over the small time
        // difference x_f = x_i + v_i * dt
        result_body->position[i] =
            original_body->position[i] +
            ( original_body->velocity[i] * TIME_DIFFERENCE );
    }

    // the result body is returned by reference
}


void outputPerformanceSummary( JobInfo* info, float total_simulation_time )
{
    // Note: putting the results to the standard out will cause them to be
    // appended to an output file by the Sun Grid Engine

    // report the data in an order that will be nice for processing later
/* Number of Nodes, Number of Bodies, Number of Simulation Steps (Time Ticks),
   Time Step Magnitude (s), Total Time to Compute Simulation (s),
   Name of File Data Was From
 */
    printf( "%d, %d, %d, %E, %E, %s\r\n",
            info->total_nodes,
            info->num_bodies,
            info->num_simulation_steps,
            info->time_difference,
            total_simulation_time,
            info->input_file_name );

    // no return - void
}


void printUniverse( JobInfo* info, Body* universe )
{
    int i = 0;
    for( i = 0; i < info->num_bodies; i++ )
    {
        printf( "Body %d\n"
                "==========\n"
                "Mass: %d\n"
                "Position:  %f | %f | %f\n"
                "Velocity:  %f | %f | %f\n",
                universe[i].position[0], universe[i].position[1], universe[i].position[2],
                universe[i].velocity[0], universe[i].velocity[1], universe[i].velocity[2] );
    }
}


