#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define G 6.67384e-11

inline F(ma, mb, r){ return  G*ma*b/r*r; }

int time;

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

int main(int* argv, char** argc)
{
    Body *BODIES = new Body[argv[0]];
    size = atoi(argv[1]);
    int i, j;

    for


    return 0;
}
