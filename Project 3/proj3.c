#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <omp.h>

#ifndef BIG
	#define BIG 1000000000
#endif

#ifndef ELEMENTS
	#define ELEMENTS 4
#endif

#ifndef THREAD_COUNT
	#define THREAD_COUNT 1
#endif 

// Set FIX to a nonzero value to enable "fixed" code and call twoFix
#ifndef FIX
	#define FIX 0
#endif

#ifndef PAD_COUNT
	#define PAD_COUNT 0
#endif

/* A struct of arbitrary values and some padding ints */
struct s
{
	float value;
	int pad[ PAD_COUNT ];
};

/* Function Prototypes */
float oneFix( struct s * );
float twoFix( struct s * );

int main(int argc, char *argv[])
{
	// Allocate array with an address at an address that is a multiple of 64
	struct s *array = aligned_alloc( (size_t)64, ( sizeof(struct s) * ELEMENTS ) );

	float (*pFun)( struct s *);	// function pointer because why not?

	// Set function pointer
	if ( FIX ) pFun = &twoFix;
	else pFun = &oneFix;

	// Print diagnostic as CSV.
	fprintf( stdout, "%d,%d,%d,%lf\n", FIX, THREAD_COUNT, PAD_COUNT, pFun( array ));

	free( array );

	return EXIT_SUCCESS;	// Way to go, (thread) team!
}

//	Performs calculations and reports MegaCalculations per second 
float oneFix( struct s *array )
{	
	omp_set_num_threads( THREAD_COUNT );
	double tStart = omp_get_wtime();

	#pragma omp parallel for
	for( int i = 0; i < ELEMENTS; i++ )
	{
		for( int j = 0; j < BIG; j++ )
		{
			// False sharing ahoy!
			array[ i ].value = array[ i ].value + 2.;
		}
	}

	double tEnd = omp_get_wtime();

	return 4000. / ( tEnd - tStart );	// MCalcs/sec
}

/*	Performs calculations and reports MegaCalculations per second
	Uses a private variable on each thread's stack. */
float twoFix( struct s *array )
{
	omp_set_num_threads( THREAD_COUNT );
	double tStart = omp_get_wtime( );

	#pragma omp parallel for
	for( int i = 0; i < ELEMENTS; i++ )
	{
		// Private variable goes on each thread's own stack.
		float tmp = array[ i ].value;

		for( int j = 0; j < BIG; j++ )
		{
			tmp = tmp + 2.;
		}

		array[ i ].value = tmp;	// Record result.
	}

	double tEnd = omp_get_wtime();
	return 4000. / ( tEnd - tStart );	// MCalcs/sec
}
