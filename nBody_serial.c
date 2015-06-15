// gcc -O1 -o nBody_serial nBody_serial.c -lm -lrt
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <math.h>

typedef double data_t;

#define TIMESTEP 			100
#define NUMBER_OF_PARTICLES	5
#define MINVAL  			0.0
#define MAXVAL  			1000.0
#define MIN_DISTANCE 		0.0001	
#define MASS 				3000        
#define GIG 				1000000000
#define CPG 				2.53           	// Cycles per GHz -- Adjust to your computer
#define EPSILON MIN_DISTANCE * MIN_DISTANCE;

typedef struct 							//structure of a particle
{
	data_t m;
	data_t x, y, z, vx, vy, vz;
} particle_struct, *particle_ptr;

typedef struct 							//structure to hold temp values of a particle
{
	data_t xnew, ynew, znew;
} temp_struct, *temp_ptr;

int initializeParticleArray(particle_ptr part_arr, long int len);
void nBody_serial_implementation(particle_ptr part_arr, temp_ptr temp_arr, long int nParticles);

main(int argc, char *argv[])
{
	long int nParticles = 0;
	int ret = 0;

	if (argc > 1) 
	{
		nParticles  = atoi(argv[1]);
	}
	else 
	{
		nParticles = NUMBER_OF_PARTICLES;
	}

	particle_ptr particleArray = (particle_ptr) malloc(nParticles * sizeof(particle_struct));
	temp_ptr tempArray = (temp_ptr) malloc(nParticles * sizeof(temp_struct));

	ret = initializeParticleArray(particleArray, nParticles);

	if(ret == 0)
		nBody_serial_implementation(particleArray, tempArray, nParticles);


	return 0;
}

/* Define fRand for random floating point */
double fRand(double fMin, double fMax)
{
    double f = (double)random()/ (double)(RAND_MAX);
    return fMin + f * (fMax - fMin);
}


int initializeParticleArray(particle_ptr part_arr, long int len)
{
	long int i = 0;

	if(len > 0)
	{
		for(i = 0; i < len; i++)
		{
			part_arr[i].x = (data_t)(fRand((double)(MINVAL),(double)(MAXVAL)));
			part_arr[i].y = (data_t)(fRand((double)(MINVAL),(double)(MAXVAL)));
			part_arr[i].z = (data_t)(fRand((double)(MINVAL),(double)(MAXVAL)));

			part_arr[i].vx = 0.0;
			part_arr[i].vy = 0.0;
			part_arr[i].vz = 0.0;

			part_arr[i].m = (data_t) MASS;
		}

		return 0;
	}

	else
		return -1;
}

void nBody_serial_implementation(particle_ptr part_arr, temp_ptr temp_arr, long int nParticles)
{
	long int i, j, k;				//iterators
	
	data_t acc_x, acc_y, acc_z;	// acceleration accumulators in x, y and z direction
	data_t dx, dy, dz;			// distance vectors in x, y and z directions

	data_t dist, distCubed;		// r factor in the gravitational force equation

	data_t force;

	data_t dt = (data_t) TIMESTEP;
	data_t dt2 = dt * dt;
	data_t r = 0;

	for(k = 0; k < 10; k++)
	{				
		printf("Iteration %ld\n", k);

		for(i = 0; i < nParticles; i++)
		{
			acc_x = 0.0;
			acc_y = 0.0;
			acc_z = 0.0;

			for(j = 0; j < nParticles; j++)
			{
				dx = part_arr[j].x - part_arr[i].x;
				dy = part_arr[j].y - part_arr[i].y;
				dz = part_arr[j].z - part_arr[i].z;

				r = dx*dx + dy*dy + dz*dz + EPSILON;

				dist = (data_t) 1.0 / sqrt(r);
				
				distCubed = dist * dist * dist;

				force = part_arr[j].m * distCubed;

				acc_x += force * dx;
				acc_y += force * dy;
				acc_z += force * dz;
			}

			temp_arr[i].xnew = part_arr[i].x + dt * part_arr[i].vx + 0.5 * dt2 * acc_x;	// New position of particle
			temp_arr[i].ynew = part_arr[i].y + dt * part_arr[i].vy + 0.5 * dt2 * acc_y;
			temp_arr[i].znew = part_arr[i].z + dt * part_arr[i].vz + 0.5 * dt2 * acc_z;

			part_arr[i].vx += dt * acc_x;
			part_arr[i].vy += dt * acc_y;
			part_arr[i].vz += dt * acc_z;	
		}

		for(i = 0; i < nParticles; i++)				// update the particle array
		{
			part_arr[i].x = temp_arr[i].xnew;
			part_arr[i].y = temp_arr[i].ynew;
			part_arr[i].z = temp_arr[i].znew;
		}

		for(i = 0; i < nParticles; i++)
		{
			printf("Particle %ld at %f,%f,%f moving with speed %f,%f,%f\n", i, part_arr[i].x, part_arr[i].y, part_arr[i].z, part_arr[i].vx, part_arr[i].vy, part_arr[i].vz);
		}
	}


}