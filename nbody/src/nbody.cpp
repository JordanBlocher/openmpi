#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define G 6.67384e-11

int TIME, NBODIES;

struct Body
{
    float mass;
    float position[3];
    float velocity[3];    
};

struct Force
{
    float x;
    float y;
};

void ComputeForces(int i)
{
    int j, r, r2, rx, ry, reciprocalForce;
    Force f;
    for(j=0; j<NBODIES; ++j)
    {
        if( i== j) continue;
        rx = universe[i].x-universe[j].x;
        ry = universe[i].y-universe[j].y;
        r2 = sqr(rx) + sqr(ry);
        r  = sqrt((double) r2);
        reciprocalForce= G * universe[i].mass * universe[j].mass
            / distance2;
        f.magnitude[0] += ReciprocalForce * rx/r;
        f.magnitude[1] += ReciprocalForce * ry/r;
        f.magnitude[2] += ReciprocalForce * rz/r;
    }
}

int main(int* argv, char** argc)
{
    TIME = atoi(argv[1]);
    int i, j, t;

    Body *temp = new Body[NBODIES];
    for(t=0; t<TIME; ++t)
    {
        for(i=0; i<NBODIES; ++i)
        {
            // Force routine
            for(j=0; j<3; ++j)
            {
                ComputeForce(i);
                temp[i].position[j] = universe[i].position[j] + Forces[i];   
            }
        } 
        delete [] universe;
        universe = temp;
    }


    return 0;
}
