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
#define DT 100 // Big time step so to see movement

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
    int type = 0;
    float mass;
    float position[3]; // center of mass
    Cell cell;
};

// Create new tree branch helper function
void initNodes(Node &universe)
{
    universe.type = 1;
    universe.mass = 0;
    universe.cell.nodes[0] = NULL;
    universe.cell.nodes[1] = NULL;
    universe.cell.nodes[2] = NULL;
    universe.cell.nodes[3] = NULL;
    universe.cell.nodes[4] = NULL;
    universe.cell.nodes[5] = NULL;
    universe.cell.nodes[6] = NULL;
    universe.cell.nodes[7] = NULL;

}

void printBodies(Body *&bodies, int t)
{
    int i;
    for(i=0; i<NBODIES; ++i)
        printf("time %d, node %p: %f %f %f %d\n", t, &bodies[i], bodies[i].position[0], bodies[i].position[1], bodies[i].position[2], bodies[i].mass);

}

void printNode(Node &node)
{
    int i;
    printf("node %p type %d: %f %f %f %f\n", &node, node.type, node.position[0], node.position[1], node.position[2], node.mass);
    if(node.type == 1)
    {
        for(i=0; i<8; i++)
        {
            if(node.cell.nodes[i] != NULL)
                printf("node[%d] type %d at %p: %f %f %f %f\n", i, node.cell.nodes[i]->type, node.cell.nodes[i], node.cell.nodes[i]->position[0], node.cell.nodes[i]->position[1], node.cell.nodes[i]->position[2], node.cell.nodes[i]->mass);
            else
                printf("node[%d] at %p: \n", i, &node.cell.nodes[i]);
        }
    }
}

void printBody(Body &body)
{
    printf("body %p: %f %f %f %d\n", &body, body.position[0], body.position[1], body.position[2], body.mass);
}

void printTree(Node *root, int level)
{
    int i;
    printf("level %d at %p: %f %f %f %f\n", level, &root, root->position[0], root->position[1], root->position[2], root->mass);
    for(i=0; i<8; i++)
    {
        if(root->cell.nodes[i] != NULL)
        {
            printf("nodes[%d] type %d at %p: %f %f %f %f\n", i, root->cell.nodes[i]->type, &(root->cell.nodes[i]), root->cell.nodes[i]->position[0], root->cell.nodes[i]->position[1], root->cell.nodes[i]->position[2], root->cell.nodes[i]->mass);
            if(root->cell.nodes[i]->type == 1)
                printTree(root->cell.nodes[i], ++level);
        }
        else
            printf("nodes[%d] type %d at %p:\n", i, root->type, &root->cell.nodes[i]); 
    }

}

void deallocateTree(Node *node) 
{
    int i;

	if ((node == NULL) || (node->type == 0))
		return;
	else 
	{
	    for(i=0; i<8; i++)
		    deallocateTree(node->cell.nodes[i]);

		free(node);
	}
}

// Read initial data and allocate memory
void readNbodyData(char* file_name, Node*& universe, Body*& bodies, Force*& forces)
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
    universe = (Node*)calloc( NBODIES, sizeof( Node ) );
    forces = (Force*)calloc( NBODIES, sizeof( Force ) );
    bodies = (Body*)calloc( NBODIES, sizeof( Body ) );
    for( i = 0; i < NBODIES; i++ )
        memcpy(forces[i].magnitude, zeros, 3*sizeof(float));

    // read in all body information
    for( i = 0; i < NBODIES; ++i )
    {
        // Temp body to allow same pointer in tree and array
        body = (Body*)calloc( NBODIES, sizeof( Body ) );

        universe[i].type = 0;
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
        universe[i].type = 0;
        universe[i].cell.body = &bodies[i];
        universe[i].mass = bodies[i].mass;
        for(j=0; j<3; j++)
            universe[i].position[j] = bodies[i].position[j];
    }

    free(body);
    // return the universe array by reference
    // (note: we are not error checking for production efficiency)
}

// Leapfrog routine
void Leapfrog(Node &body, Force &force)
{
    int j,k;
    Body* temp = (Body*)calloc(1, sizeof( Body ));
    float vminushalf[3], vplushalf[3];

        for(j=0; j<3; ++j)
        {
            for(k=0; k<3; ++k)
            {
                // Leapfrog : v(t - 1/2)
                vminushalf[k] = body.cell.body->velocity[k];   
            }
            // Leapfrog : v(t + 1/2)
            vplushalf[j] = 0.5*(body.cell.body->velocity[j] + force.magnitude[j]/body.mass*DT);
            // v(t)
            temp->velocity[j] = (vplushalf[j] - vminushalf[j])*body.cell.body->mass*DT;   
            // x(t + 1/2)
            temp->position[j] = body.cell.body->position[j] + vplushalf[j]*DT;   
            temp->mass = body.cell.body->mass;
        }
    body.cell.body = temp;

    free(temp);
}

// Force routine
void ComputeForce(Node &mi, Node &mj, Force &force) 
{
    int k;
    float d, r2, r[3], reciprocalForce;
    // Compute actions of body i on the universe
    r2 = 0;
    for(k=0; k<3; k++)
    {
        r[k] = mi.position[k] - mj.position[k]; //rx, ry, rz
        r2 += r[k]*r[k]; // r*r
    }
    d  = sqrt((double) r2); // distance = sqrt(rx2 + ry2 + rz2)
    // F = G*mj*mi/r*r (gravitational force between body i & j)
    reciprocalForce = G * mi.mass * mj.mass/ r2;
    for(k=0; k<3; k++)
        force.magnitude[k] = reciprocalForce * r[k]/d; // Action of body i on body j
}

// Recursive tree force computation
void ComputeForceRecursive(Node &system, Node &body, Force &force, float d2)
{
    int j;
    float r[3], r2=0;

    for(j=0; j<3; j++)
    {
        r[j] = system.position[j] - body.position[j]; //rx, ry, rz
        r2 += r[j]*r[j]; // r*r
    }
    if(r2 < d2) // Recurse through smaller solar systems
    {
        if(system.type == 1)
        {
            d2*=0.25; // Decrease diameter of system for cutoff approximation
            for(j=0; j<8; ++j)
            {
                if(system.cell.nodes[j] != NULL)
                    ComputeForceRecursive(&(*system.cell.nodes[j]), body, force, d2);
                else continue;
            }
        }
        else // Compute force for local body
        {
            ComputeForce(system, body, force);
        }
    }
    else // Compute force approximation for solar system as mass cluster
    {
        ComputeForce(system, body, force);
    }
}

// Compute center of mass of tree or branch
void ComputeCOM(Node &node) 
{
    int m = 0;
    float tempPosition[3] = {0, 0, 0};

    Node *temp;

    int j = 0;
    node.mass = m; // Initialize to zero
    int i;
    for (i = 0; i < 8; i++) 
    {
        temp = node.cell.nodes[i];
        if (temp != NULL) 
        {
            // Swap out nodes to sum mass (going down the tree)
            node.cell.nodes[i] = NULL;
            node.cell.nodes[j] = temp;
            j++;

            if (temp.type != 0) // Recurse if not at a leaf
                ComputeCOM(temp);

            m = temp.mass;
            node.mass += m;
            tempPosition[0] += temp.position[0] * m;
            tempPosition[1] += temp.position[1] * m;
            tempPosition[2] += temp.position[2] * m;
        }
    }

    // Place this node at the center of the system
    m = 1.0 / node.mass;
    for(i=0; i<3; i++)
        node.position[i] = tempPosition[i] * m;
    
    free(temp);
}

// Compute diameter of node and center
void ComputeDTC(Node &universe, float *center, float &diameter) // compute distance to center
{
    float min[3] = { 1.0e6, 1.0e6, 1.0e6 };
    float max[3] = { -1.0e6, -1.0e6, -1.0e6 }; 
    float position[3];

    int i, j;
    // Find largest distance between bodies
    for (i = 0; i < NBODIES; i++) 
    {
        position[0] = universe[i].position[0];
        position[1] = universe[i].position[1];
        position[2] = universe[i].position[2];

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

void insert(Node &system, Node &node, float r) 
{
	bool finished = false;
    int j, idx1, idx2;
    float rprev;
    Node *temp, *newSystem;
	do 
    {
		float d[3];
        idx1=0;
        for(j=0; j<3; j++)
        {
            if (system.position[j] < node.position[j]) 
            {
                // Part of tree algorithm, why?? A bit confused here, compute the index where the new node goes 
                idx1 = (j == 0) ? 1 : (j == 1) ? idx1+2 : idx1+4;
                d[j] = r; // Set maximum distance to be in system
            }
            else d[j] = 0; // If new system zero is max distance
        }
        // We have room, make the node a child 
		if (system.cell.nodes[idx1] == NULL) 
        {
			system.cell.nodes[idx1] = node;
			finished = true;
		} 
        // System exists here, don't insert node, go further into the tree
        else if (system.cell.nodes[idx1].type == 1)
        {
			r *= 0.5; // Make new system smaller
			system = system.cell.nodes[idx1];
		} 
        // Tree is not deep enough, make the current system a branch 
        else 
        {
			rprev = 0.5 * r;
			newSystem = (Node*)calloc(1, sizeof(Node)); // New node as branch root
            initNodes(newSystem);
            for(j=0; j<3; j++)
                newSystem.position[j] = system.position[j] - rprev + d[j]; 
            idx2 = 0;
            for(j=0; j<3; j++)
			    if (newSystem.position[j] < node.position[j]) 
                    idx2 = (j == 0) ? 1 : (j == 1) ? idx2+2 : idx2+4;
            
            // Swap out crowded node with branch
			newSystem.cell.nodes[idx2] = node;
			temp = system.cell.nodes[idx1];
			system.cell.nodes[idx1] = newSystem;
			system = newSystem;
			node = temp;
			r = rprev;
		}
	} while (!finished);
}

// Create Body type for message-passing
void CreateMPIDatatype()
{
    // Array of datatypes
    MPI_Datatype type[5] = { MPI_LB, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_UB };
    // Number of each type
    int block_len[5] = {1, 1, 3, 3, 1};
    
    // Displacements
    MPI_Aint disp[5];
    Body body[2]; // <-- used to calculate displacement to next body

    // MPI gets addresses for displacement in memory buffer
    MPI_Get_address(&body[0], &disp[0]);
    MPI_Get_address(&(body[0].mass), &disp[1]);
    MPI_Get_address(&(body[0].position), &disp[2]);
    MPI_Get_address(&(body[0].velocity), &disp[3]);
    MPI_Get_address(&(body[1].velocity), &disp[4]);

    disp[1] = disp[1] - disp[0];
    disp[2] = disp[2] - disp[0];
    disp[3] = disp[3] - disp[0];
    disp[4] = disp[4] - disp[0];

    // MPI Body struct
    MPI_Type_create_struct(5, block_len, disp, type, &MPI_Body);

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


void updateBodies()
{
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
    Node *universe;
    Force *forces;
    Body *bodies;
    Node *root;

    std::chrono::high_resolution_clock::time_point startclock, stopclock;
    std::chrono::duration<double> time_span;

    // Init MPI & get numprocs and rank
    MPI_Init(&argc,&argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    CHECKMPI(err); 

    CreateMPIDatatype(); 

    std::string filename = "data/init.dat";
    readNbodyData((char*)filename.c_str(), universe, bodies, forces);

    printBodies(bodies, -1);
    cout<<"\n";
    //MPI_Bcast body pointer array
    //MPI_Bcast(bodies, NBODIES, MPI_Body, 0, MPI_COMM_WORLD);
    //CHECKMPI(err); 

    for(t=0; t<TIME; ++t)
    {
        cout<<"\nAllocating root\n";
        root = (Node*)calloc(1, sizeof( Node ));
        cout<<"\n..root address " <<root<<" ok?\n";
        printBodies(bodies, t);
        cout<<"\n";
        // Compute diameter and center of universe
        ComputeDTC(universe, center, diameter);
        // Set root node
        initNodes(root);
        radius = diameter * 0.5;
  
        // Set root to center of universe  
        for(j=0; j<3; j++)
            root.position[j] = center[j];
	
        // Build Octree
        for (j = 0; j < NBODIES; ++j) 
			insert(root, universe[j], radius); 
        printTree(root, 0);
  
        // Compute center of mass of the universe
        ComputeCOM(root);

        // Force routines
        for (i = 0; i < NBODIES; ++i) 
	       ComputeForceRecursive(root, universe[i], forces[i], diameter*diameter);
        
        //MPI_Allgather body pointer array
        //MPI_Allgather(bodies, NBODIES, MPI_Body, bodies, NBODIES, MPI_Body, MPI_COMM_WORLD);

        for (i = 0; i < NBODIES; ++i) 
           Leapfrog(universe[i], forces[i]);
          
        // Delete tree
        deallocateTree(root);
    }


    if(rank == 0)
    {
        stopclock = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(stopclock - startclock);

        //printf("\n %d %d %f\n",numprocs, N, (double)time_span.count());
    }

    MPI_Finalize();
    return 0;
}
