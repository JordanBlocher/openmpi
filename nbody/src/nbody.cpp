#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstring>

using std::cout;
using std::endl;

#define TIME_DIFFERENCE 0.0003
#define G 6.67384e-11
#define DT 0.003

int TIME, NBODIES;

struct Body
{
    int index;
    float mass;
    float position[3];
    float velocity[3];    
};

struct Force
{
    float magnitude[3];  // each magnitude corresponds to a cartesian direction
};

struct Universe
{
    float COM[3]; 
};

struct Octree
{
    float COM[3]; // center of mass
    Octree *nodes[8];
    Universe *universe; 
};

void readNbodyData( char* file_name, Body**& universe );

void readNbodyData( char* file_name, Body**& universe )
{
  // variables 
  FILE* file = NULL;
  int i = 0;

  // open the file
  file = fopen( file_name );

  // read in the number of bodies
  fscanf( file, "%d", &NBODIES );

  // allocate the array to contain pointers to all of the bodies
  universe = calloc( NBODIES, sizeof( Body* ) );
  for( i = 0; i < NBODIES; i++ )
  {
    // allocate the bodies themselves
    universe[i] = calloc( 1, sizeof( Body ) );
  }

  // read in all body information
  for( i = 0; i < NBODIES; ++i )
  {
    // read in the mass value
    fscanf( file, "%f", &(universe[i]->mass) );

    // read in the body's position as a spacial coordinate triple
    fscanf( file, "%f %f %f", &(universe[i]->position[0]),
            &(universe[i]->position[1]), &(universe[i]->position[2]) );
  }

  // return the universe array by reference
  // (note: we are not error checking for production efficiency)
}


void ComputeForce(int i, Body **universe, Force **forces) // wrong way to pass universe, hold please
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
        reciprocalForce = G * universe[i]->mass * universe[j]->mass/ r2; // F = G*(ma - mb)/r*r
        for(k=0; k<3; k++)
            forces[j]->magnitude[k] += reciprocalForce * r[k]/d; // Action of body i on body j
    }
}

int main(int, char** argv)
{
    TIME = atoi(argv[1]);
    int i, j, t;
    float vminushalf[3], vplushalf[3];
    Body *universe;
    Force *forces;

    // Temp universe for pointer swapping
    Body *temp = new Body[NBODIES];
    for(t=0; t<TIME; ++t)
    {
        for(i=0; i<NBODIES; ++i)
        {
            // Force routine
            for(j=0; j<3; ++j)
            {
                // Leapfrog : v(t - 1/2)
                vminushalf[j] = universe[i].velocity[j];   
            }
            ComputeForce(i, &universe, &forces);
            for(j=0; j<3; ++j)
            {
                // Leapfrog : v(t + 1/2)
                vplushalf[j] = 0.5*(universe[i].velocity[j] + forces[i].magnitude[j]/universe[i].mass*DT);
                // v(t)
                temp[i].velocity[j] = (vplushalf[j] - vminushalf[j])*universe[i].mass*DT;   
                // x(t + 1/2)
                temp[i].position[j] = universe[i].position[j] + universe[i].position[j] + vplushalf[j]*DT;   
            }
        } 
        delete [] universe;
        universe = temp;
    }


    return 0;
}
