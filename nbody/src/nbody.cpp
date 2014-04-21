#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstring>
#include <sstream>

using std::cout;
using std::endl;

#define TIME_DIFFERENCE 0.0003
#define G 1.0 //6.67384e-2  // Factor all values by 1 billion
#define DT 1.0 // Big time step so to see movement

int TIME, NBODIES;
float zeros[3] = {0, 0, 0};

struct Body
{
    int index;
    int mass;
    float position[3];
    float velocity[3];    
};

struct Force
{
    float magnitude[3];  // each magnitude corresponds to a cartesian direction
};

struct Octree
{
    float COM[3]; // center of mass
    Octree *nodes[8];
    Body *universe; 
};

void readNbodyData( char* file_name, Body**& universe, Force**& forces)
{
    // variables 
    FILE* file = NULL;
    int i = 0;

    // open the file
    file = fopen( file_name , "r");

    // read in the number of bodies
    fscanf( file, "%d", &NBODIES );

    // allocate the array to contain pointers to all of the bodies
    universe = (Body**)calloc( NBODIES, sizeof( Body* ) );
    forces = (Force**)calloc( NBODIES, sizeof( Force* ) );
    for( i = 0; i < NBODIES; i++ )
    {
        // allocate the bodies themselves
        universe[i] = (Body*)calloc( 1, sizeof( Body ) );
        // allocate memory for forces
        forces[i] = (Force*)calloc( 1, sizeof( Force ) );
        memcpy(forces[i]->magnitude, zeros, 3*sizeof(float));
    }

    // read in all body information
    for( i = 0; i < NBODIES; ++i )
    {
        // read in the body's position as a spacial coordinate triple
        fscanf( file, "%f %f %f", &(universe[i]->position[0]),
                &(universe[i]->position[1]), &(universe[i]->position[2]) );

        // read in the initial velocity for the body
        fscanf( file, "%f %f %f", &(universe[i]->velocity[0]),
                &(universe[i]->velocity[1]), &(universe[i]->velocity[2]) );

        // read in the mass value
        fscanf( file, "%d", &(universe[i]->mass) );
        universe[i]->velocity[0] *= 0;
        universe[i]->velocity[1] *= 0;
        universe[i]->velocity[2] *= 0;
    }

    // return the universe array by reference
    // (note: we are not error checking for production efficiency)
}


void ComputeForce(int i, Body **universe, Force **&forces) 
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
            r[k] = universe[i]->position[k]-universe[j]->position[k]; //rx, ry, rz
            r2 += r[k]*r[k]; // r*r
        }
        d  = sqrt((double) r2); // distance = sqrt(rx2 + ry2 + rz2)
        // F = G*mj*mi/r*r (gravitational force between body i & j)
        reciprocalForce = G * universe[i]->mass * universe[j]->mass/ r2;
        for(k=0; k<3; k++)
            forces[j]->magnitude[k] = reciprocalForce * r[k]/d; // Action of body i on body j
    }
}

int main(int argc, char** argv)
{
    if(argc < 2) 
    {
        cout <<" Usage: nbody <time>\n\n";
        return 0;
    }
 
    TIME = atoi(argv[1]);
    int i, j, k, t;
    float vminushalf[3], vplushalf[3];
    Body **universe, **temp, **swap;
    Force **forces;
    FILE* file = NULL;
    std::ostringstream ss;

    readNbodyData("data/init.dat", universe, forces);

    // Temp universe for pointer swapping
    temp = (Body**)calloc( NBODIES, sizeof( Body* ) );
    for(i=0; i<NBODIES; ++i)
        temp[i] = (Body*)calloc(1, sizeof(Body));

    for(t=0; t<TIME; ++t)
    {
        // Print state of universe
        ss << "data/state" << t << ".dat";
        file = fopen(ss.str().c_str(), "w");
        for(i=0; i<NBODIES; ++i)
            fprintf( file, "%f %f %f %d\n", universe[i]->position[0], universe[i]->position[1], universe[i]->position[2], universe[i]->mass);
        ss.str("");


            // Force routine
        for(i=0; i<NBODIES; ++i)
            ComputeForce(i, universe, forces);
        for(i=0; i<NBODIES; ++i)
        {
             for(j=0; j<3; ++j)
            {
                for(k=0; k<3; ++k)
                {
                    // Leapfrog : v(t - 1/2)
                    vminushalf[k] = universe[i]->velocity[k];   
                }
                // Leapfrog : v(t + 1/2)
                vplushalf[j] = 0.5*(universe[i]->velocity[j] + forces[i]->magnitude[j]/universe[i]->mass*DT);
                // v(t)
                temp[i]->velocity[j] = (vplushalf[j] - vminushalf[j])*universe[i]->mass*DT;   
                //temp[i]->velocity[j] = universe[i]->velocity[j] + forces[i]->magnitude[j]*DT/(float)universe[i]->mass;   
                // x(t + 1/2)
                temp[i]->position[j] = universe[i]->position[j] + vplushalf[j]*DT;   
                //temp[i]->position[j] = universe[i]->position[j] + temp[i]->velocity[j]*DT;   
                temp[i]->mass = universe[i]->mass;
            }
        }

        swap = universe;
        universe = temp;
        temp = swap;
      }


      return 0;
                }
