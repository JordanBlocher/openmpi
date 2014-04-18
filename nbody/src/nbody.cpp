#define CHECKMPI(err) if(err){ fprintf(stderr,"Communicator failed.\n");MPI_Abort(MPI_COMM_WORLD,1);};

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

using std::cout;
using std::endl;

#define TIME_DIFFERENCE 0.0003
#define G 6.67384e-2  // Factor all values by 1 billion
#define DT 100

int TIME, NBODIES;
float zeros[3] = {0, 0, 0};
int numprocs, rank, idx;

MPI_Datatype MPI_Body;

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

// Read initial data and allocate memory
void readNbodyData( char* file_name, Node**& universe, Body**& bodies, Force**& forces)
{
    // variables 
    FILE* file = NULL;
    int i, j;
    Body *body;

    // open the file
    file = fopen( file_name , "r");

    // read in the number of bodies
    fscanf( file, "%d", &NBODIES );

    // allocate the array to contain pointers to all of the bodies
    universe = (Node**)calloc( NBODIES, sizeof( Node* ) );
    forces = (Force**)calloc( NBODIES, sizeof( Force* ) );
    bodies = (Body**)calloc( NBODIES, sizeof( Body* ) );
    body = (Body*)calloc( NBODIES, sizeof( Body ) );
    for( i = 0; i < NBODIES; i++ )
    {
        // allocate the bodies themselves
        universe[i] = (Node*)calloc( 1, sizeof( Node ) );
        // allocate memory for forces
        forces[i] = (Force*)calloc( 1, sizeof( Force ) );
        memcpy(forces[i]->magnitude, zeros, 3*sizeof(float));
    }

    // read in all body information
    for( i = 0; i < NBODIES; ++i )
    {
        universe[i]->type = 0;
        // read in the body's position as a spacial coordinate triple
        fscanf( file, "%f %f %f", &(body->position[0]),
                &(body->position[1]), &(body->position[2]) );

        // read in the initial velocity for the body
        fscanf( file, "%f %f %f", &(body->velocity[0]),
                &(body->velocity[1]), &(body->velocity[2]) );

        // read in the mass value
        fscanf( file, "%d", &(body->mass) );

        bodies[i] = body;

        universe[i]->type = 0;
        universe[i]->cell.body = bodies[i];
        universe[i]->mass = bodies[i]->mass;
        for(j=0; j<3; j++)
            universe[i]->position[j] = bodies[i]->position[j];

        free(body);

    }

    // return the universe array by reference
    // (note: we are not error checking for production efficiency)
}

// Force routine
void ComputeForce(int i, Node**& universe, Force**& forces) 
{
    int j, k;
    float d, r2, r[3], reciprocalForce;
    // Compute actions of body i on the universe
    for(j=0; j<NBODIES; ++j)
    {
        if( i== j) continue;
        r2 = 0;
        for(k=0; k<3; k++)
        {
            r[k] = universe[i]->cell.body->position[k] - universe[j]->cell.body->position[k]; //rx, ry, rz
            r2 += r[k]*r[k]; // r*r
        }
        d  = sqrt((double) r2); // distance = sqrt(rx2 + ry2 + rz2)
        // F = G*mj*mi/r*r (gravitational force between body i & j)
        reciprocalForce = G * universe[i]->mass * universe[j]->mass/ r2;
        for(k=0; k<3; k++)
            forces[j]->magnitude[k] = reciprocalForce * r[k]/d; // Action of body i on body j
    }
}

// Recursive tree force computation
void ComputeForceRecursive(Node**& universe, Force**& forces, Body**& bodies, int t)
{
    int i, j; 
    float vminushalf[3], vplushalf[3];
    for(i=0; i<NBODIES; ++i)
        ComputeForce(i, universe, forces);
    for(i=0; i<NBODIES; ++i)
    {
        for(j=0; j<3; ++j)
        {
            for(j=0; j<3; ++j)
            {
                // Leapfrog : v(t - 1/2)
                vminushalf[j] = universe[i]->velocity[j];   
            }
            // Leapfrog : v(t + 1/2)
            vplushalf[j] = 0.5*(universe[i]->velocity[j] + forces[i]->magnitude[j]/universe[i]->mass*DT);
            // v(t)
            temp[i]->velocity[j] = (vplushalf[j] - vminushalf[j])*universe[i]->mass*DT;   
            // x(t + 1/2)
            temp[i]->position[j] = universe[i]->position[j] + universe[i]->position[j] * vplushalf[j]*DT;   
            temp[i]->mass = universe[i]->mass;
        }
    }
}

// Compute center of mass
void ComputeCOM(Node *&node) // compute center of mass
{
    int m = 0, p[3];

    Node *ch;

    int j = 0;
    node->mass = m;
    int i;
    for (i = 0; i < 8; i++) 
    {
        ch = node->cell.nodes[i];
        if (ch != NULL) 
        {
            node->cell.nodes[i] = NULL;
            node->cell.nodes[j] = ch;
            j++;

            if (ch->type != 0)
                ComputeCOM(ch);

            m = ch->mass;
            node->mass += m;
            p[0] += ch->position[0] * m;
            p[1] += ch->position[1] * m;
            p[2] += ch->position[2] * m;
        }
    }

    m = 1.0 / node->mass;
    node->position[0] = zeros[0] * m;
    node->position[1] = zeros[1] * m;
    node->position[2] = zeros[2] * m;
}

// Compute diameter of node and center
void ComputeDTC(Node**& universe, float *center, float diameter) // compute distance to center
{
    float min[3] = { 1.0e6, 1.0e6, 1.0e6 };
    float max[3] = { -1.0e6, -1.0e6, -1.0e6 }; 
    float position[3];

    int i, j;
    for (i = 0; i < NBODIES; i++) 
    {
        position[0] = universe[i]->position[0];
        position[1] = universe[i]->position[1];
        position[2] = universe[i]->position[2];

        for(j=0; j< 3; ++j)
        {
            if (min[j] > position[j])
                min[j] = position[j];

            if (max[j] < position[j])
                max[j] = position[j];
        }

    }
    diameter = max[0] - min[0];

    for(j=1; j< 3; ++j)
    {
        if (diameter < (max[j] - min[j]))
            diameter = (max[j] - min[j]);
    }

    for(j=0; j< 3; ++j)
        center[j] = (max[j] + min[j]) * 0.5;

}

void CreateMPIDatatype()
{
    MPI_Datatype type[5] = { MPI_LB, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_UB };
    int block_len[5] = {1, 1, 3, 3, 1};
    MPI_Aint disp[5];
    Body body[2];

    MPI_Get_address(&body[0], &disp[0]);
    MPI_Get_address(&(body[0].mass), &disp[1]);
    MPI_Get_address(&(body[0].position), &disp[2]);
    MPI_Get_address(&(body[0].velocity), &disp[3]);
    MPI_Get_address(&(body[1].velocity), &disp[4]);

    disp[1] = disp[1] - disp[0];
    disp[2] = disp[2] - disp[0];
    disp[3] = disp[3] - disp[0];
    disp[4] = disp[4] - disp[0];

    MPI_Type_create_struct(5, block_len, disp, type, &MPI_Body);

    MPI_Type_commit(&MPI_Body);
}

// Print state of universe
void printBodiesToFile(Body**& bodies, int t)
{
    int i;
    FILE* file = NULL;
    std::ostringstream ss;
    ss << "data/state" << t << ".dat";
    file = fopen(ss.str().c_str(), "w");
    for(i=0; i<NBODIES; ++i)
        fprintf( file, "%f %f %f %d\n", bodies[i]->position[0], bodies[i]->position[1], bodies[i]->position[2], bodies[i]->mass);
    ss.str("");

}

void initBodies(Node*& universe)
{
    universe->type = 1;
    universe->mass = 0;
    universe->cell.nodes[0] = NULL;
    universe->cell.nodes[1] = NULL;
    universe->cell.nodes[2] = NULL;
    universe->cell.nodes[3] = NULL;
    universe->cell.nodes[4] = NULL;
    universe->cell.nodes[5] = NULL;
    universe->cell.nodes[6] = NULL;
    universe->cell.nodes[7] = NULL;

}

int main(int argc, char** argv)
{
    if(argc < 2) 
    {
        cout <<" Usage: nbody <time>\n\n";
        return 0;
    }

    TIME = atoi(argv[1]);
    int i, j, t, err;
    float center[3], diameter=0, radius=0;;
    Node **universe, **temp, **swap;
    Force **forces;
    Body **bodies;
    idx = 0;

    std::chrono::high_resolution_clock::time_point startclock, stopclock;
    std::chrono::duration<double> time_span;

    std::string filename = "data/init.dat";
    readNbodyData((char*)filename.c_str(), universe, bodies, forces);

    // Temp universe for pointer swapping
    temp = (Node**)calloc( NBODIES, sizeof( Node* ) );
    for(i=0; i<NBODIES; ++i)
        temp[i] = (Node*)calloc(1, sizeof(Node));

    // Init MPI & get numprocs and rank
    //MPI_Init(&argc,&argv);
    //err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //CHECKMPI(err); 

    for(t=0; t<TIME; ++t)
    {
        initBodies(*universe);
        ComputeDTC(universe, center, diameter);
        radius = diameter * 0.5;
        for(j=0; j<3; j++)
            (*universe)->position[i] = center[i];

		for (j = 0; j < NBODIES; ++j) 
		{
			insert(universe, bodies[j], radius); 
		}

		ComputeCOM(*universe);

		for (j = 0; j < NBODIES; ++j) 
		{
			ComputeForceRecursive(*universe, *forces[j], *bodies[j], diameter, t);
		}
        swap = universe;
        universe = temp;
        temp = swap;
    }


    if(rank == 0)
    {
        stopclock = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stopclock - startclock);

        //printf("\n %d %d %f\n",numprocs, N, (double)time_span.count());
    }

    //MPI_Finalize();
    return 0;
}
