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
#define DT 10 // Big time step so to see movement

int TIME, NBODIES;
float zeros[3] = {0, 0, 0};
int numprocs, rank;

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
    int mass;
    float position[3]; // center of mass
    Cell cell;
};

// Create new tree branch helper function
void initNodes(Node *universe)
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

void printBody(Body *body)
{
    printf("body %p: %f %f %f %d\n", body, body->position[0], body->position[1], body->position[2], body->mass);
}

void printBodies(Body *&bodies, int)
{
    int i;
    for(i=0; i<NBODIES; ++i)
    {
        printBody(&bodies[i]);
    }

}

void printUniverse(Node **&universe, int)
{
    int i;
    for(i=0; i<NBODIES; ++i)
    {
        printBody(universe[i]->cell.body);
    }

}

void printTree(Node *root, int level)
{
    int i;
    printf("level %d type %d at %p: %f %f %f %d\n", level, root->type, root, root->position[0], root->position[1], root->position[2], root->mass);
    if(root->type == 0)
        printBody(root->cell.body);
    else if(root->type == 1)
    { 
        for(i=0; i<8; i++)
        {
            if(root->cell.nodes[i] != NULL)
            {
                printf("nodes[%d] type %d at %p: %f %f %f %d\n", i, root->cell.nodes[i]->type, &(root->cell.nodes[i]), root->cell.nodes[i]->position[0], root->cell.nodes[i]->position[1], root->cell.nodes[i]->position[2], root->cell.nodes[i]->mass);
                if(root->cell.nodes[i]->type == 1)
                    printTree(&(*root->cell.nodes[i]), ++level);
            }
        }
    }

}

// Read initial data and allocate memory
void readNbodyData(char* file_name, Node**& universe, Body*& bodies, Force*& forces)
{
    // variables 
    FILE* file = NULL;
    int i, j;
    Body* body;

    // open the file
    file = fopen( file_name , "r");

    // read in the number of bodies
    fscanf( file, "%d", &NBODIES );

    // allocate the array to contain pointers to all of the bodies
    universe = (Node**)calloc( NBODIES, sizeof( Node* ) );
    forces = (Force*)calloc( NBODIES, sizeof( Force ) );
    bodies = (Body*)calloc( NBODIES, sizeof( Body ) );

    // read in all body information
    for( i = 0; i < NBODIES; ++i )
    {
        universe[i] = (Node*)calloc( NBODIES, sizeof( Node ) );
        memcpy(forces[i].magnitude, zeros, 3*sizeof(float));

        // Temp body to allow same pointer in tree and array
        body = (Body*)calloc(1, sizeof( Body ) );

        universe[i]->type = 0;
        // read in the body's position as a spacial coordinate triple
        fscanf( file, "%f %f %f", &(body->position[0]),
                &(body->position[1]), &(body->position[2]) );

        // read in the initial velocity for the body
        fscanf( file, "%f %f %f", &(body->velocity[0]),
                &(body->velocity[1]), &(body->velocity[2]) );

        // read in the mass value
        fscanf( file, "%d", &(body->mass) );

        bodies[i] = *body; 

        // Body pointer also stored in universe for binary space recursion
        universe[i]->type = 0;
        universe[i]->cell.body = &bodies[i];
        universe[i]->mass = universe[i]->cell.body->mass;
        for(j=0; j<3; j++)
            universe[i]->position[j] = universe[i]->cell.body->position[j];
    }

    // return the universe array by reference
    // (note: we are not error checking for production efficiency)
}

// Leapfrog routine
void Leapfrog(Node *body, Force *force)
{
    int j,k;
    Body* temp = body->cell.body;
    float vminushalf[3], vplushalf[3];

    for(j=0; j<3; ++j)
    {
        for(k=0; k<3; ++k)
        {
            // Leapfrog : v(t - 1/2)
            vminushalf[k] = body->cell.body->velocity[k];   
        }
        // Leapfrog : v(t + 1/2)
        vplushalf[j] = 0.5*(body->cell.body->velocity[j] + force->magnitude[j]/body->mass*DT);
        // v(t)
        temp->velocity[j] = (vplushalf[j] - vminushalf[j])*body->cell.body->mass*DT;   
        // x(t + 1/2)
        temp->position[j] = body->cell.body->position[j] + vplushalf[j]*DT;   
        temp->mass = body->cell.body->mass;
    }

}

// Force routine
void F(Node *mi, Node *mj, Force *force)
{
    int k;
    float d, r2, r[3], reciprocalForce;

    // Compute actions of body i on the universe
    r2 = 0;
    for(k=0; k<3; k++)
    {
        r[k] = mi->position[k] - mj->position[k]; //rx, ry, rz
        r2 += r[k]*r[k]; // r*r
    }

    d  = sqrt((double) r2); // distance = sqrt(rx2 + ry2 + rz2)
    // F = G*mj*mi/r*r (gravitational force between body i & j)
    reciprocalForce = G * mi->mass * mj->mass/ r2;

    for(k=0; k<3; k++)
        force->magnitude[k] += reciprocalForce * r[k]/d; // Action of body i on body j
}

// Recursive tree force computation
void ComputeForceRecursive(Node *system, Node *body, Force *force, float d2)
{
    int j;
    float r[3], r2=0;

    for(j=0; j<3; j++)
    {
        r[j] = system->position[j] - body->position[j]; //rx, ry, rz
        r2 += r[j]*r[j]; // r*r
    }
    if(r2 < d2) // Recurse through smaller solar systems
    {
        if(system->type == 1)
        {
            d2*=0.25; // Decrease diameter of system for cutoff approximation
            for(j=0; j<8; ++j)
            {
                if(system->cell.nodes[j] != NULL)
                    ComputeForceRecursive(&(*system->cell.nodes[j]), &(*body), &(*force), d2);
                else continue;
            }
        }
        else if (system != body)// Compute force for local body
        {
            F(&(*system), &(*body), &(*force));
        }
    }
    else // Compute force approximation for solar system as mass cluster
    {
        F(&(*system), &(*body), &(*force));
    }
}

// Recursive Sum
void ComputeForce(Node *system, Node *body, Force *force, float d2) 
{
    int k;

    for(k=0; k<3; k++)
        force->magnitude[k] = 0;

    ComputeForceRecursive(&(*system), &(*body), &(*force), d2);
}


// Compute center of mass of tree or branch
void ComputeCOM(Node *node) 
{
    int m = 0;
    float tempPosition[3] = {0, 0, 0};

    Node *temp;

    int j = 0;
    node->mass = m; // Initialize to zero
    int i;
    for (i = 0; i < 8; i++) 
    {
        temp = node->cell.nodes[i];
        if (temp != NULL) 
        {
            // Swap out nodes to sum mass (going down the tree)
            node->cell.nodes[i] = NULL;
            node->cell.nodes[j] = temp;
            j++;

            if (temp->type != 0) // Recurse if not at a leaf
                ComputeCOM(temp);

            m = temp->mass;
            node->mass += m;
            tempPosition[0] += temp->position[0] * m;
            tempPosition[1] += temp->position[1] * m;
            tempPosition[2] += temp->position[2] * m;
        }
    }

    // Place this node at the center of the system
    m = 1.0 / node->mass;
    for(i=0; i<3; i++)
        node->position[i] = tempPosition[i] * m;
}

// Compute diameter of node and center
void ComputeDTC(Node**& universe, float *center, float &diameter) // compute distance to center
{
    float min[3] = { 1.0e6, 1.0e6, 1.0e6 };
    float max[3] = { -1.0e6, -1.0e6, -1.0e6 }; 
    float position[3];

    int i, j;
    // Find largest distance between bodies
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
    diameter = max[0] - min[0]; // Default diameter to x-axis

    // Compute largest diameter in x, y, z
    for(j=1; j< 3; ++j)
    {
        if (diameter < (max[j] - min[j]))
            diameter = (max[j] - min[j]);
    }

    // The center is the midpoint
    for(j=0; j< 3; ++j)
        center[j] = (max[j] + min[j]) * 0.5;

}

void insert(Node *system, Node *node, float r) 
{
    int j, idx1, idx2;
    float rprev;
    Node *newSystem, *temp;
    float d[3] = {0, 0, 0};
    idx1=0;
    //printf("System %p \n", system);
    //printf("Node %p pos %f %f %f\n", node, node->position[0],node->position[1],node->position[2]);
    for(j=0; j<3; j++)
    {
        if (system->position[j] < node->position[j]) 
        {
            // Part of tree algorithm, why?? A bit confused here, compute the index where the new node goes 
            idx1 = (j == 0) ? 1 : (j == 1) ? idx1+2 : idx1+4;
            d[j] = r; // Set maximum distance to be in system
        }
    }
    // We have room, make the node a child 
    if (system->cell.nodes[idx1] == NULL) 
    {
        system->cell.nodes[idx1] = node;
        return;
    } 
    // System exists here, don't insert node, go further into the tree
    else if (system->cell.nodes[idx1]->type == 1)
    {
        //printf("%p type %d idx %d\n", system->cell.nodes[idx1],system->cell.nodes[idx1]->type,idx1);
        insert(&(*system->cell.nodes[idx1]), &(*node), r * 0.5);
    } 
    // Tree is not deep enough, make the current system a branch 
    else 
    {
        rprev = 0.5 * r;
        newSystem = (Node*)calloc(1, sizeof(Node)); // New node as branch root
        initNodes(newSystem);
        for(j=0; j<3; j++)
            newSystem->position[j] = system->position[j] - rprev + d[j]; 
        //printf("oldSys pos %f %f %f\n", system->position[0],system->position[1],system->position[2]);
        //printf("newSys pos %f %f %f\n", newSystem->position[0],newSystem->position[1],newSystem->position[2]);
        //printf("node pos %f %f %f\n", node->position[0],node->position[1],node->position[2]);
        idx2 = 0;
        for(j=0; j<3; j++)
            if (newSystem->position[j] < node->position[j]) 
                idx2 = (j == 0) ? 1 : (j == 1) ? idx2+2 : idx2+4;

        // Swap out crowded node with branch
        newSystem->cell.nodes[idx2] = node;
        temp = &(*system->cell.nodes[idx1]);
        system->cell.nodes[idx1] = newSystem;
        node = temp;
        //printf("new system %p\n", newSystem);
        //printf("Branched node %p \n", system->cell.nodes[idx1]);
        insert(&(*newSystem), &(*temp), rprev);
    }
}

// Create Body type for message-passing
void CreateMPIDatatype()
{
    // Array of datatypes
    MPI_Datatype type[3] = {  MPI_FLOAT, MPI_FLOAT, MPI_FLOAT };
    // Number of each type
    int block_len[3] = { 1, 3, 3 };
    // Displacements
    MPI_Aint disp[3] = { 0, sizeof(float),
        ( 4 * sizeof(float) ) };
    // MPI Body struct
    MPI_Type_create_struct(3, block_len, disp, type, &MPI_Body);
    MPI_Type_commit(&MPI_Body);
}

// Print state of universe
void printBodiesToFile(Body*& bodies, int t)
{
    int i;
    FILE* file = NULL;
    std::ostringstream ss;
    ss << "data/state" << t << ".dat";
    file = fopen(ss.str().c_str(), "w");
    for(i=0; i<NBODIES; ++i)
        fprintf( file, "%f %f %f %d\n", bodies[i].position[0], bodies[i].position[1], bodies[i].position[2], bodies[i].mass);
    ss.str("");

}

int main(int argc, char** argv)
{
    if(argc < 3) 
    {
        cout <<" Usage: nbody <time> <filename>\n\n";
        return 0;
    }

    TIME = atoi(argv[1]);
    int i, j, k, t, err, nlocal_bodies;
    float center[3], diameter=0, radius=0;;
    Node **universe;
    Force *forces;
    Body *bodies;
    Body *local_bodies;
    Node *root;

    std::chrono::high_resolution_clock::time_point startclock, stopclock;
    std::chrono::duration<double> time_span;

    // Init MPI & get numprocs and rank
    MPI_Init(&argc,&argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    CHECKMPI(err); 

    CreateMPIDatatype(); 
    std::string filename = argv[2];
    readNbodyData((char*)filename.c_str(), universe, bodies, forces);

    nlocal_bodies = NBODIES / numprocs; 

    if(rank == 0)
    {
        //nlocal_bodies = NBODIES / numprocs; 

        // MPI_Bcast body data
        //err = MPI_Bcast(&NBODIES, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //CHECKMPI(err); 

        //err = MPI_Bcast(bodies, NBODIES, MPI_Body, 0, MPI_COMM_WORLD);
        //CHECKMPI(err); 

        for(t=0; t<TIME; ++t)
        {
            cout<<"TIME "<<t<<endl;
            root = (Node*)calloc(1, sizeof( Node ));
            cout<<"Bodies\n";
            printBodies(bodies, t);
            cout<<"\n\n";

            // Compute diameter and center of universe
            ComputeDTC(universe, center, diameter);
            // Set root node
            initNodes(root);
            radius = diameter * 0.5;

            // Set root to center of universe  
            for(j=0; j<3; j++)
                root->position[j] = center[j];

            // Build Octree
            for (j = 0; j < NBODIES; ++j) 
                insert(root, universe[j], radius); 

            // Compute center of mass of the universe
            ComputeCOM(&(*root));

            // Force routines
            for (i = 0; i < nlocal_bodies; ++i) 
                ComputeForce(&(*root), &(*universe[i]), &(forces[i]), diameter*diameter);

            //MPI_Allgather body pointer array
            err = MPI_Allgather(bodies, nlocal_bodies, MPI_Body, bodies, nlocal_bodies, MPI_Body, MPI_COMM_WORLD);
            CHECKMPI(err);

            for (i = 0; i < NBODIES; ++i) 
                Leapfrog(&(*universe[i]), &(forces[i]));
        }

        stopclock = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stopclock - startclock);

        printf("\n %d %d %f\n", numprocs, NBODIES, (double)time_span.count());
    }

    else if (rank != 0)
    {
        // Get Body data
        //err = MPI_Bcast(&NBODIES, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //CHECKMPI(err); 
        
        //nlocal_bodies = NBODIES / numprocs; 
        //universe = (Node**)calloc( NBODIES, sizeof( Node* ) );
        //for(i=0; i<NBODIES; i++)
         //   universe[i] = (Node*)calloc(1, sizeof( Node ));
        //forces = (Force*)calloc( NBODIES, sizeof( Force ) );
        //bodies = (Body*)calloc( NBODIES, sizeof( Body ) );

        //err = MPI_Bcast(bodies, NBODIES, MPI_Body, 0, MPI_COMM_WORLD);
        //CHECKMPI(err); 

        int start = (rank * nlocal_bodies);
        int stop = start + nlocal_bodies;

        for(t=0; t<TIME; ++t)
        {
            root = (Node*)calloc(1, sizeof( Node ));

            // Compute diameter and center of universe
            ComputeDTC(universe, center, diameter);
            // Set root node
            initNodes(root);
            radius = diameter * 0.5;

            // Set root to center of universe  
            for(j=0; j<3; j++)
                root->position[j] = center[j];

            // Build Octree
            for (j = 0; j < NBODIES; ++j) 
            {
                //universe[j]->type = 0;
                //universe[j]->cell.body = &bodies[j];
                //universe[j]->mass = bodies[j].mass;
                //for(k=0; k<3; k++)
                  //  universe[j]->position[k] = bodies[j].position[k];
                insert(&(*root), &(*universe[j]), radius); 
            }

            // Compute center of mass of the universe
            ComputeCOM(&(*root));
            // Force routines
            for (i = start; i < stop; ++i) 
                ComputeForce(&(*root), &(*universe[i]), &(forces[i]), diameter*diameter);

            //MPI_Allgather body pointer array
            local_bodies = &bodies[start];
            err = MPI_Allgather(local_bodies, nlocal_bodies, MPI_Body, bodies, nlocal_bodies, MPI_Body, MPI_COMM_WORLD);
            CHECKMPI(err);

            for (i = 0; i < NBODIES; ++i) 
                Leapfrog(&(*universe[i]), &(forces[i]));
        }
    }

    MPI_Finalize();
    return 0;
}
