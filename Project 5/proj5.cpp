
#include "simd.p5.h"
#include <cfloat>
#include <functional>

using std::function; 	// GetTime parameters
using std::fill_n; 		// Array filling

// Set to nonzero for verbose logging
#ifndef VERBOSE
	#define VERBOSE 0
#endif

/* float array element count
	Set to 0 to print header and exit */
#ifndef LEN
	#define LEN 1000
#endif

/* Number of loops to average to ensure data integrity
	Set to > 1 to enable confirmation */
#ifndef TRIES
	#define TRIES 1
#endif

// Number of attempts to get good data before aborting
#ifndef ABORT
	#define ABORT 256
#endif

// Global arrays.
float A[LEN];
float B[LEN]; 
float C[LEN];

// Function prototypes
void err( const char * );
void SisdMul( float *, float *, float *, int );
float SisdMulSum( float *a, float *b, float *c, int len );
float GetTime( function<void( float *, float *, float *, int )> );
float GetTime( function<float( float *, float *, int )> );

int main (int argc, char **argv)
{
	// Fill arrays with arbitrary values
	fill_n( B, LEN, 222.323 );
	fill_n( C, LEN, 323.222 );

	// Get the data
	float 	SpeedupMul = 	GetTime( SisdMul ) 		/ GetTime( SimdMul ),
			SpeedupMulSum =	GetTime( SisdMulSum )	/ GetTime( SimdMulSum );

	// Display the header and exit on zero LEN
	if ( !( LEN ) ) {
		fprintf ( stdout, "Array Size,SIMD Mult Speedup,SIMD Mult+Red Speedup\n" );
		return EXIT_SUCCESS;
	}

	// Display the data
	fprintf( stdout, "%d,%8.2lf,%8.2lf\n", LEN, SpeedupMul, SpeedupMulSum );

	return EXIT_SUCCESS;
}

// Print error message and exit with error
void err( const char *msg )
{
	fprintf( stderr, "%s\n", msg );
	exit( EXIT_FAILURE );
}

// Single Instruction Single Data array multiplication
void SisdMul( float *a, float *b, float *c, int len )
{
	for( int i = 0; i < len; i++ ) {
		c[i] = a[i] * b[i];
	}
}

// Single Instruction Single Data array multiplication + reduction
float SisdMulSum( float *a, float *b, float *c, int len )
{
	float sum = 0.;

	for( int i = 0; i < len; i++ ) {
		sum += a[i] * b[i];
	}

	return sum;
}

// Gets a verified reliable time for peak array multiply performance.
float GetTime( function<void(float *, float *, float *, int)> Mul )
{
	int 	count = 0;

	double	time0,
			time1,
			ratio,
		 	avg = 0,
			peak = DBL_MAX;

	do{
		count++;

		for (int i = 0; i < TRIES; i++ ) {
			time0 = omp_get_wtime();
			Mul( A, B, C, LEN );
			time1 = omp_get_wtime();

			// Record times
			avg += time1 - time0;
			if ( peak > time1 - time0 ) peak = time1 - time0;
		}

		avg /= TRIES;

		// Calculate ratio for peak & avg
		ratio = peak > avg ? avg / peak : peak / avg;
	
	// Repeat until peak is within 20% of average
	} while ( (ratio <= .8) && count < ABORT );

	if (count >= ABORT) err("No reliable data available. Aborted!");

	if ( VERBOSE ) fprintf( stderr, "Count: %d\tavg: %lf\tpeak: %lf\tratio: %lf\n", count, avg, peak, ratio);

	return peak;
}

/* Gets a verified reliable time for peak array multiply + reduce performance.
	Yes, there is duplicate code. I didn't take the time to handle differing signatures. :-| */
float GetTime( function<float(float *, float *, int)> MulSum )
{
	int 	count = 0;

	double	time0,
			time1,
			ratio,
		 	avg = 0,
			peak = DBL_MAX;

	do{
		count++;

		for (int i = 0; i < TRIES; i++ ) {
			time0 = omp_get_wtime();
			MulSum( A, B, LEN ); // We don't care about the return value
			time1 = omp_get_wtime();

			// Record times
			avg += time1 - time0;
			if ( peak > time1 - time0 ) peak = time1 - time0;
		}

		avg /= TRIES;

		// Calculate ratio for peak & avg
		ratio = peak > avg ? avg / peak : peak / avg;

		
	// Repeat until peak is within 20% of average
	} while ( (ratio <= .8) && count < ABORT );

	if (count >= ABORT) err("No reliable data available. Aborted!");

	if ( VERBOSE ) fprintf ( stderr, "Count: %d\tavg: %lf\tpeak: %lf\tratio: %lf\n", count, avg, peak, ratio);

	return peak;
}