#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

// Default NUMT to serial if not compiled with -D
#ifndef NUMT
	#define NUMT		1
#endif

// Default NUMNODES to 512 if not compiled with -D
#ifndef NUMNODES 
	#define NUMNODES 512
#endif

// Default NUMTRIES to 8 if not compiled with -D
#ifndef NUMTRIES
	#define NUMTRIES
#endif

// Surface definition
#define XMIN	 0.
#define XMAX	 3.
#define YMIN	 0.
#define YMAX	 3.

#define TOPZ00  0.
#define TOPZ10  1.
#define TOPZ20  0.
#define TOPZ30  0.

#define TOPZ01  1.
#define TOPZ11  6.
#define TOPZ21  1.
#define TOPZ31  0.

#define TOPZ02  0.
#define TOPZ12  1.
#define TOPZ22  0.
#define TOPZ32  4.

#define TOPZ03  3.
#define TOPZ13  2.
#define TOPZ23  3.
#define TOPZ33  3.

#define BOTZ00  0.
#define BOTZ10  -3.
#define BOTZ20  0.
#define BOTZ30  0.

#define BOTZ01  -2.
#define BOTZ11  10.
#define BOTZ21  -2.
#define BOTZ31  0.

#define BOTZ02  0.
#define BOTZ12  -5.
#define BOTZ22  0.
#define BOTZ32  -6.

#define BOTZ03  -3.
#define BOTZ13   2.
#define BOTZ23  -8.
#define BOTZ33  -3.

// Surface height function
double Height( int iu, int iv )	// iu,iv = 0 .. NUMNODES-1
{
	double u = (double)iu / (double)(NUMNODES-1);
	double v = (double)iv / (double)(NUMNODES-1);

	// the basis functions:

	double bu0 = (1.-u) * (1.-u) * (1.-u);
	double bu1 = 3. * u * (1.-u) * (1.-u);
	double bu2 = 3. * u * u * (1.-u);
	double bu3 = u * u * u;

	double bv0 = (1.-v) * (1.-v) * (1.-v);
	double bv1 = 3. * v * (1.-v) * (1.-v);
	double bv2 = 3. * v * v * (1.-v);
	double bv3 = v * v * v;

	// finally, we get to compute something:

    double top =		bu0 * ( bv0*TOPZ00 + bv1*TOPZ01 + bv2*TOPZ02 + bv3*TOPZ03 )
					+ bu1 * ( bv0*TOPZ10 + bv1*TOPZ11 + bv2*TOPZ12 + bv3*TOPZ13 )
					+ bu2 * ( bv0*TOPZ20 + bv1*TOPZ21 + bv2*TOPZ22 + bv3*TOPZ23 )
					+ bu3 * ( bv0*TOPZ30 + bv1*TOPZ31 + bv2*TOPZ32 + bv3*TOPZ33 );

	double bot =		bu0 * ( bv0*BOTZ00 + bv1*BOTZ01 + bv2*BOTZ02 + bv3*BOTZ03 )
					+ bu1 * ( bv0*BOTZ10 + bv1*BOTZ11 + bv2*BOTZ12 + bv3*BOTZ13 )
					+ bu2 * ( bv0*BOTZ20 + bv1*BOTZ21 + bv2*BOTZ22 + bv3*BOTZ23 )
					+ bu3 * ( bv0*BOTZ30 + bv1*BOTZ31 + bv2*BOTZ32 + bv3*BOTZ33 );

	return top - bot;	// if the bottom surface sticks out above the top surface
						// then that contribution to the overall volume is negative
}


int main() 
{

#ifndef _OPENMP
	fprintf( stderr, "OpenMP is not supported here -- sorry.\n" );
	return 1;
#endif

	double	volume,			// Volume result
			mhpsAvg,			// Megaheights per second
			mhpsPeak,			// Megaheights per second
			tStart,			// Wall time start
			tEnd,			// Wall time end
			tLast,			// Duration of most recent calculation
			tLastAvg,
			tLastPeak,
			tSerial,		// Duration of single threaded loop
			tSerialAvg,
			tSerialPeak,
			speedupAvg,		
			speedupPeak,		
			efficiencyAvg,
			efficiencyPeak,
			fParallelAvg,
			fParallelPeak;
		
	// Print csv headers	
	fprintf( stdout, "Nodes,Threads,Volume,Avg Time,Avg Mh/s,Avg Speedup,Avg Efficiency,Avg Fp,");
	fprintf( stdout, "Peak Time,Peak Mh/s,Peak Speedup,Peak Efficiency Avg,Peak Fp\n");

	// Test NUMNODES against threadcounts from 1 to NUMT
	for( int curT = 1; curT <= NUMT; curT++ )
	{

		omp_set_num_threads( curT );

		tLastAvg = 0;
		tLastPeak = DBL_MAX; 

		for ( int i = 0; i < NUMTRIES; i++)

		{

		volume = 0;

		tStart = omp_get_wtime( );

		#pragma omp parallel for default(none),reduction(+:volume)
		for( int i = 0; i < NUMNODES*NUMNODES; i++ )
		{
			// Get current coordinates
			int iu = i % NUMNODES;
			int iv = i / NUMNODES;
		
			// Check for edge/corner & get area coefficient: edge = 0.5, corner = 0.25
			double edgeFactor = 1.;
			if ( iu == 0 || iu == ( NUMNODES -1 ) ) edgeFactor *= 0.5;
			if ( iv == 0 || iv == ( NUMNODES -1 ) ) edgeFactor *= 0.5;
	
			// Calculate area if tile is full-sized	
			double fullTileArea =	(  ( ( XMAX - XMIN )/(double)(NUMNODES-1) )  *
								( ( YMAX - YMIN )/(double)(NUMNODES-1) )  );

			// Calculate actual volume of column
			volume += edgeFactor * fullTileArea * Height( iu, iv );

		}

		tEnd = omp_get_wtime( );

		// Calculate duration
		tLast = tEnd - tStart;

		// Accumulate for avg
		tLastAvg += tLast;

		// Record peaks
		if ( tLast < tLastPeak ) tLastPeak = tLast;

		}

		// Calculate average
		tLastAvg /= (double) NUMTRIES;

		// Calculate Megaheights per second
		mhpsAvg = pow(NUMNODES, 2.) / (tLastAvg) / 1000000.;
		mhpsPeak = pow(NUMNODES, 2.) / (tLastPeak) / 1000000.;

		// Record serial performance
		if ( curT == 1 ) 
		{
			tSerialAvg = tLastAvg; 
			tSerialPeak = tLastPeak; 
			speedupAvg = 0;
			efficiencyAvg = 0;
			fParallelAvg = 0;
			speedupPeak = 0;
			efficiencyPeak = 0;
			fParallelPeak = 0;
		}
		// Record parallel performance
		else
		{
			speedupAvg = tSerialAvg / tLastAvg;
			efficiencyAvg = speedupAvg / curT;
			fParallelAvg = (curT / (curT - 1)) * ((tSerialAvg - tLastAvg) / tSerialAvg );

			speedupPeak = tSerialPeak / tLastPeak;
			efficiencyPeak = speedupPeak / curT;
			fParallelPeak = (curT / (curT - 1)) * ((tSerialPeak - tLastPeak) / tSerialPeak );
		}

		fprintf( stdout, "%d,%d,%lf,", NUMNODES, curT, volume );
		fprintf( stdout, "%lf,%lf,%lf,%lf,%lf,", tLastAvg, mhpsAvg, speedupAvg, efficiencyAvg, fParallelAvg );
		fprintf( stdout, "%lf,%lf,%lf,%lf,%lf\n", tLastPeak, mhpsPeak, speedupPeak, efficiencyPeak, fParallelPeak );

	}

	return 0;	

}


