/*
	Name: Joshua L. Rogers
	Course: CS475 Sp2018
	Assignment: Project 4
	Date: 14 May 2018

	Description:	Month-by-month simulation of a grain-growing operation.
					The amount the grain grows is affected by the temperature,
					amount of precipitation, and the number of graindeer
					around to eat it. The number of graindeer depends on the
					amount of grain available to eat.
*/

#include <stdlib.h>	// for EXIT_FAILURE / EXIT_SUCCESS
#include <stdio.h> 	// for printing results
#include <math.h>
#include <omp.h>	// OpenMP
#include <time.h>	// time() for seeding pRNG
#include <stdbool.h> // bool type

#ifndef M_PI
	#define M_PI 3.14159265358979323846264338327950288
#endif

// Option macros to be set on compile if desired.

// Minimum 3 required for basic funtionality
// Set to 4 to enable student-added section
#ifndef THREADCOUNT
	#define THREADCOUNT 3
#endif

// Simulation begins this month (Default: January)
#ifndef MONTHSTART
	#define MONTHSTART 0
#endif

// Simulation begins this year
#ifndef YEARSTART
	#define YEARSTART 2017
#endif

// Simulation stops when this year is reached.
#ifndef YEARSTOP
	#define YEARSTOP 2023
#endif

// Initial grain height
#ifndef GRAINSTART
	#define GRAINSTART 1.
#endif

/* Initial deer population
		Fun Fact: Graindeer reproduce asexually via spores.
*/ 
#ifndef DEERSTART
	#define DEERSTART 1
#endif

// Initial gas over norm 
#ifndef GASSTART
	#define GASSTART 0
#endif

/*
	Global State Variables
		length units: inches
		temperature units: degF
*/
int	NowMonth = MONTHSTART;		// 0 - 11
int	NowYear = YEARSTART;		// 2017 - 2022
int	NowNumDeer = DEERSTART;		// number of deer in the current population
float NowHeight = GRAINSTART;	// grain height in inches
float NowGreenGas = GASSTART;	// additional atmospheric greenouse gas content

float NowPrecip;	// inches of rain per month
float NowTemp;		// temperature this month

const float GRAIN_GROWS_PER_MONTH = 8.0;
const float ONE_DEER_EATS_PER_MONTH = 0.5;

const float AVG_PRECIP_PER_MONTH = 6.0;	// average
const float AMP_PRECIP_PER_MONTH = 6.0;	// variance / amplitude
const float RANDOM_PRECIP = 2.0;	// plus or minus noise
const float MIDPRECIP = 10.0;

const float AVG_TEMP = 50.0;	// average
const float AMP_TEMP = 20.0;	// variance / amplitude
const float RANDOM_TEMP = 10.0;	// plus or minus noise
const float MIDTEMP = 40.0; 


void GrainDeer();
void Grain();
void Watcher( unsigned int *seedp );
void Greenhouse();
void SetWeather();
void PrintState( bool state );
float Ranf( unsigned int *seedp,  float low, float high );
float SQR( float x );
float in2cm ( float x );
float F2C ( float x );


int main (int argc, char **argv)
{
	// Check for invalid starting state
	if (THREADCOUNT < 3 || NowYear >= YEARSTOP || NowMonth < 0 || NowMonth > 11) exit(EXIT_FAILURE);

	// Set weather for initial month
	unsigned int seed = time(NULL);
	SetWeather( &seed );

	// Print header
	PrintState( false );

	omp_set_num_threads( THREADCOUNT );
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			GrainDeer( );
		}

		#pragma omp section
		{
			Grain( );
		}

		#pragma omp section
		{
			Watcher( &seed );
		}

		#if THREADCOUNT > 3
			#pragma omp section
			{
				Greenhouse( );
			}
		#endif
	}
}

// Calculates graindeer population
void GrainDeer()
{
	while(NowYear < YEARSTOP) {

		// Compute new graindeer count in local variable
		int NewNumDeer = NowNumDeer;
		if ( (float) NowNumDeer < NowHeight) NewNumDeer++;
		else if ( (float) NowNumDeer > NowHeight) NewNumDeer--;
		#pragma omp barrier // Done Computing

		// Copy new graindeer count into global variable
		NowNumDeer = NewNumDeer;
		#pragma omp barrier	// Done Assigning
		#pragma omp barrier // Done Printing
	}
}

// Calculates grain height
void Grain()
{
	while(NowYear < YEARSTOP) {

		// Calculate divergence from ideal weather
		float tempFactor = exp(   -SQR(  ( NowTemp - MIDTEMP ) / 10.  )   );
		float precipFactor = exp(   -SQR(  ( NowPrecip - MIDPRECIP ) / 10.  )   );	

		// Calculate new height into a local variable
		float NewHeight = NowHeight;
		NewHeight += tempFactor * precipFactor * GRAIN_GROWS_PER_MONTH;
		NewHeight -= (float)NowNumDeer * ONE_DEER_EATS_PER_MONTH;
		if (NewHeight < 0. ) NewHeight = 0.; // No negative heights!
		#pragma omp barrier // Done Computing

		NowHeight = NewHeight;
		#pragma omp barrier	// Done Assigning
		#pragma omp barrier // Done Printing
	}
}

// Prints results, increments month/year, sets weather
void Watcher(unsigned int *seedp)
{
	
	while(NowYear < YEARSTOP) {
		#pragma omp barrier // Done Computing
		#pragma omp barrier	// Done Assigning
		
		// PRINT RESULTS
		PrintState( true );

		// Increment Month/Year
		if ( ++NowMonth == 12 ) {
			NowMonth = 0;
			NowYear++;
		}

		// Set new weather conditions
		SetWeather(seedp);
		#pragma omp barrier // Done Printing
	}
}

/*	Computes the impact of graindeer & grain on
	atmospheric greenhouse gas content as percentage
	over/under baseline */
void Greenhouse()
{
	while(NowYear < YEARSTOP) {

		float NewGreenGas = (float) NowNumDeer; // Deer produce gas (e.g., CO2)		
		NewGreenGas -= (float) NowHeight / 2.5; // Grain absorbs gas (i.e., CO2)
		#pragma omp barrier // Done Computing

		NowGreenGas += NewGreenGas; // Update global gas value.
		if ( NowGreenGas <= -100. ) NowGreenGas = -99; // Assume not ALL greenhouse gas can be removed.

		#pragma omp barrier	// Done Assigning
		#pragma omp barrier // Done Printing
	}
}

// Set the NowTemp and NowPrecip global floats
void SetWeather(unsigned int *seedp)
{
	// Calculate impact of gas content on weather
	float GasFactor = 1. + ( NowGreenGas / 100. );

	float ang = (  30.*(float)NowMonth + 15.  ) * ( M_PI / 180. );

	 // Greenhouse effect on temperature
	float temp = ( AVG_TEMP - AMP_TEMP * cos( ang ) ) * GasFactor;
	// Greenhouse effect increases unpredictability
	NowTemp = temp + Ranf( seedp, ( -RANDOM_TEMP * GasFactor ), ( RANDOM_TEMP * GasFactor ) ); 

	// Greenhouse effect on baseline precipitation
	float precip = ( AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin( ang ) ) / GasFactor; 
	// Greenhouse effect increases unpredictability
	NowPrecip = precip + Ranf( seedp,  ( -RANDOM_PRECIP * GasFactor ), ( RANDOM_PRECIP * GasFactor ) );  

	// "Unrain" is forbidden!
	if( NowPrecip < 0. ) NowPrecip = 0.;
}

void PrintState( bool state )
{
	// If state not indicated, print header and return
	if ( !state ) {
		fprintf( stdout, "Month,Year,Temp (C),Precip (cm),Height (cm),NumDeer");

		#if THREADCOUNT > 3
			fprintf( stdout, ",CO2 Over Baseline (\%)" );
		#endif

		fprintf( stdout, "\n" );

		return;
	}

	// Otherwise print state values
	fprintf( stdout, "%d,%d,%lf,%lf,%lf,%d", NowMonth,NowYear,F2C(NowTemp),in2cm(NowPrecip),in2cm(NowHeight),NowNumDeer);

	#if THREADCOUNT > 3
		fprintf( stdout, ",%lf", NowGreenGas );
	#endif

	fprintf( stdout, "\n" );
}

// Get a random float between ilow and ihigh
float Ranf( unsigned int *seedp,  float low, float high )
{
        float r = (float) rand_r( seedp );              // 0 - RAND_MAX

        return(   low  +  r * ( high - low ) / (float)RAND_MAX   );
}
/*
// Get a random int between ilow and ihigh
int Rani( unsigned int *seedp, int ilow, int ihigh )
{
        float low = (float)ilow;
        float high = (float)ihigh + 0.9999f;

        return (int)(  Ranf(seedp, low,high) );
}
*/
// Squares a float
float SQR( float x )
{
	return x*x;
}

// Converts inches to centimeters
float in2cm ( float x )
{
	return x * 2.54;
}

// Converts degrees Farenheit to degrees Celsius
float F2C ( float x )
{
	return (5./9.)*(x-32);
}