#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TIME_DIFFERENCE 0.0003
#define G 6.67384e-11

inline Fx(ma, mb, r){ return  G*ma*b/r*r; }

int time;

struct Body
{
    float mass;
    float position[3];
    float velocity[3];    
};

struct Force
{
    float magnitude[3];  // each magnitude corresponds to a cartesian direction
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
            &(universe[i]->position[1]), &(universe[i]-position[2]) );
  }

  // return the universe array by reference
  // (note: we are not error checking for production efficiency)
}



int main(int* argv, char** argc)
{
    Body *BODIES = new Body[argv[0]];
    size = atoi(argv[1]);
    int i, j;

    for


    return 0;
}
