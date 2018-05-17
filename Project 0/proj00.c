#include <omp.h>
#include <stdio.h>
#include <math.h>

// Pass values using -D, but provide defaults:

#ifndef NUMT
	#define NUMT		1
#endif

#ifndef ARRAYSIZE
	#define ARRAYSIZE 100000 
#endif

#ifndef NUMTRIES
	#define NUMTRIES 100000
#endif

float A[ARRAYSIZE];
float B[ARRAYSIZE];
float C[ARRAYSIZE];

int main() 
{

#ifndef _OPENMP
	fprintf( stderr, "OpenMP is not supported here -- sorry.\n" );
	return 1;
#endif

	omp_set_num_threads( NUMT );
	fprintf( stderr, "Using %d threads\n", NUMT );

	double maxMegaMults = 0.;
	double sumMegaMults = 0.;
	double start = omp_get_wtime( );

	for( int t = 0; t < NUMTRIES; t++ )
	{
		double time0 = omp_get_wtime( );

		#pragma omp parallel for
		for( int i = 0; i < ARRAYSIZE; i++ )
		{
			C[i] = A[i] * B[i];
		}

		double time1 = omp_get_wtime( );
		double megaMults = (double)ARRAYSIZE/(time1-time0)/1000000.;
		sumMegaMults += megaMults;
		if( megaMults > maxMegaMults )
			maxMegaMults = megaMults;
	}

	double end = omp_get_wtime( );
	double avgMegaMults = sumMegaMults/(double)NUMTRIES;

    printf( "Elapsed time: %8.2lf Sec\n", end - start);
	printf( "   Peak Performance = %8.2lf MegaMults/Sec\n", maxMegaMults );
	printf( "Average Performance = %8.2lf MegaMults/Sec\n", avgMegaMults );

	// note: %lf stands for "long float", which is how printf prints a "double"
	//	%d stands for "decimal integer", not "double"
	
	return 0;
}
