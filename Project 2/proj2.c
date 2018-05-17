#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>

/*	Compile using -D to set NUMTHREADS, GRAIN in [0:1], OMP_SCHED in [1:4]
	Default 1 thread, coarse-grained parallelism (0), static scheduling (1)
	May also set NUMBODIES and NUMSTEPS, but not required for this exercise.
*/

#ifndef ITERATIONS
	#define ITERATIONS 1
#endif

#ifndef NUMBODIES
	#define NUMBODIES 100
#endif

#ifndef NUMSTEPS
	#define NUMSTEPS  200
#endif

#ifndef NUMTHREADS
  #define NUMTHREADS 1
#endif

/*	OMP_SCHED reflects OMP_SCHEDULE environment variable:
	1: static	(default)
	2: dynamic 
	3: guided
	4: auto	
*/

#ifndef OMP_SCHED
	#define OMP_SCHED 1
#endif

// Now set SCHEDULE, default is static

#if (OMP_SCHED == 1)
	#define SCHEDULE static
#elif (OMP_SCHED == 2)
	#define SCHEDULE dynamic
#elif (OMP_SCHED == 3)
	#define SCHEDULE guided
#elif (OMP_SCHED == 4)
	#define SCHEDULE auto
#else
	#define SCHEDULE static
#endif

/*	GRAIN reflects the granularity of parallelism:
	0: coarse-grained		(default)
	1: fine-grained
*/

#ifndef GRAIN
	#define GRAIN 0
#endif

#if GRAIN != 0
	#if GRAIN != 1
		#define GRAIN 0
	#endif
#endif

#ifndef HEADER
	#define HEADER 0
#endif

// constants:
const double G = 6.67300e-11; // m^3 / ( kg s^2 )
const double EARTH_MASS = 5.9742e24; // kg
const double EARTH_DIAMETER = 12756000.32; // meters
const double TIMESTEP = 1.0; // secs

struct body
{
  float mass;
  float x, y, z; // position
  float vx, vy, vz; // velocity
	float fx, fy, fz; // forces
	float xnew, ynew, znew;
	float vxnew, vynew, vznew;
};

typedef struct body Body;

Body Bodies[NUMBODIES];

// function prototypes:
float GetDistanceSquared( Body *, Body * );
float GetUnitVector( Body *, Body *, float *, float *, float * );
float Ranf( float, float );
int Ranf( int, int );

int main( int argc, char *argv[ ] )
{
	#ifndef _OPENMP
		fprintf( stderr, "OpenMP is not available\n" );
		return 1;
	#endif
 
	omp_set_num_threads( NUMTHREADS );

	for( int i = 0; i < NUMBODIES; i++ )
	{
		Bodies[i].mass = EARTH_MASS  * Ranf( 0.5f, 10.f );
		Bodies[i].x = EARTH_DIAMETER * Ranf( -100.f, 100.f );
		Bodies[i].y = EARTH_DIAMETER * Ranf( -100.f, 100.f );
		Bodies[i].z = EARTH_DIAMETER * Ranf( -100.f, 100.f );
		Bodies[i].vx = Ranf( -100.f, 100.f );;
		Bodies[i].vy = Ranf( -100.f, 100.f );;
		Bodies[i].vz = Ranf( -100.f, 100.f );;
	};

	double tAvg = 0;		// Holds average time for all iterations
	double tMin = DBL_MAX;	// Holds peak minimum time across all iterations

	for (int i = 0; i < ITERATIONS; i++)
	{

		double time0 = omp_get_wtime( );

		for( int t = 0; t < NUMSTEPS; t++ )
		{
			#if GRAIN == 0
				#pragma omp parallel for default(none) shared(Bodies) schedule(SCHEDULE)
			#endif
			for( int i = 0; i < NUMBODIES; i++ )
			{
				float fx = 0.;
				float fy = 0.;
				float fz = 0.;
				Body *bi = &Bodies[i];
		
				#if GRAIN == 1
					#pragma omp parallel for default(none) shared(Bodies, i, bi) reduction(+:fx,fy,fz) schedule(SCHEDULE)
				#endif			
				for( int j = 0; j < NUMBODIES; j++ )
				{
					if( j == i ) continue;

					Body *bj = &Bodies[j];

					float rsqd = GetDistanceSquared( bi, bj );

					if( rsqd > 0. )
					{
						float f = G * bi->mass * bj->mass / rsqd;
						float ux, uy, uz;
						GetUnitVector( bi, bj, &ux, &uy, &uz );
						fx += f * ux;
						fy += f * uy;
						fz += f * uz;
					}
				}

				float ax = fx / Bodies[i].mass;
				float ay = fy / Bodies[i].mass;
				float az = fz / Bodies[i].mass;

				Bodies[i].xnew = Bodies[i].x + Bodies[i].vx*TIMESTEP + 0.5*ax*TIMESTEP*TIMESTEP;
				Bodies[i].ynew = Bodies[i].y + Bodies[i].vy*TIMESTEP + 0.5*ay*TIMESTEP*TIMESTEP;
				Bodies[i].znew = Bodies[i].z + Bodies[i].vz*TIMESTEP + 0.5*az*TIMESTEP*TIMESTEP;

				Bodies[i].vxnew = Bodies[i].vx + ax*TIMESTEP;
				Bodies[i].vynew = Bodies[i].vy + ay*TIMESTEP;
				Bodies[i].vznew = Bodies[i].vz + az*TIMESTEP;
			}

			// setup the state for the next animation step:
		 
			for( int i = 0; i < NUMBODIES; i++ )
			{
			  Bodies[i].x = Bodies[i].xnew;
			  Bodies[i].y = Bodies[i].ynew;
			  Bodies[i].z = Bodies[i].znew;
			  Bodies[i].vx = Bodies[i].vxnew;
			  Bodies[i].vy = Bodies[i].vynew;
			  Bodies[i].vz = Bodies[i].vznew;
			}

		}  // t

		double time1 = omp_get_wtime( );

		tAvg += time1 - time0;

		if (time1 - time0 < tMin) tMin = time1 - time0;

	}

	tAvg /= ITERATIONS;

	float MbpsAvg = (float) (NUMBODIES * NUMBODIES * NUMSTEPS) / tAvg / 1000000.;
	float MbpsPeak = (float) (NUMBODIES * NUMBODIES * NUMSTEPS) / tMin / 1000000.;

	// Print results as csv.
	fprintf( stdout, "%d,%lf,%lf", NUMTHREADS, MbpsAvg, MbpsPeak);

	#if OMP_SCHED == 1
		fprintf( stdout, ",static" );
	#elif OMP_SCHED == 2
		fprintf( stdout, ",dynamic" );
	#elif OMP_SCHED == 3
		fprintf( stdout, ",guided" );
	#elif OMP_SCHED == 4
		fprintf( stdout, ",auto" );
	#else
		fprintf( stdout, ",static" );
	#endif

	#if GRAIN == 1
		fprintf( stdout, ",fine");
	#else
		fprintf( stdout, ",coarse");
	#endif

	fprintf( stdout, "\n");
	
	return 0;

}


float GetDistanceSquared( Body *bi, Body *bj )
{
	float dx = bi->x - bj->x;
	float dy = bi->y - bj->y;
	float dz = bi->z - bj->z;
	return dx*dx + dy*dy + dz*dz;
}

float GetUnitVector( Body *from, Body *to, float *ux, float *uy, float *uz )
{
	float dx = to->x - from->x;
	float dy = to->y - from->y;
	float dz = to->z - from->z;
	float d = sqrt( dx*dx + dy*dy + dz*dz );

	if( d > 0. )
	{
		dx /= d;
		dy /= d;
		dz /= d;
	}
 
	*ux = dx;
	*uy = dy;
	*uz = dz;
 
	return d;
}

float Ranf( float low, float high )
{
	float r = (float) rand();         // 0 - RAND_MAX
	return(   low  +  r * ( high - low ) / (float)RAND_MAX   );
}

int Ranf( int ilow, int ihigh )
{
	float low = (float)ilow;
	float high = (float)ihigh + 0.9999f;
	return (int)(  Ranf(low,high) );
}
