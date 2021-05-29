#include <omp.h>
#include <xmmintrin.h>
#include <stdio.h>
#include <math.h>
#include <values.h>
#define SSE_WIDTH		4

#ifndef ARRAYSIZE
#define ARRAYSIZE	1000
#endif
#ifndef NUMT
#define NUMT	1
#endif
// omp_set_num_threads( NUMT );	// set the number of threads to use in parallelizing the for-loop:`
// time0 = omp_get_wtime( );
void SimdMul( float *a, float *b,   float *c,   int len ){
	int limit = ( len/SSE_WIDTH ) * SSE_WIDTH;
	register float *pa = a;
	register float *pb = b;
	register float *pc = c;
	for( int i = 0; i < limit; i += SSE_WIDTH )
	{
		_mm_storeu_ps( pc,  _mm_mul_ps( _mm_loadu_ps( pa ), _mm_loadu_ps( pb ) ) );
		pa += SSE_WIDTH;
		pb += SSE_WIDTH;
		pc += SSE_WIDTH;
	}

	for( int i = limit; i < len; i++ )
	{
		c[i] = a[i] * b[i];
	}
}


float SimdMulSum( float *a, float *b, int len ) {
	float sum[4] = { 0., 0., 0., 0. };
	int limit = ( len/SSE_WIDTH ) * SSE_WIDTH;
	register float *pa = a;
	register float *pb = b;

	__m128 ss = _mm_loadu_ps( &sum[0] );
	for( int i = 0; i < limit; i += SSE_WIDTH )
	{
		ss = _mm_add_ps( ss, _mm_mul_ps( _mm_loadu_ps( pa ), _mm_loadu_ps( pb ) ) );
		pa += SSE_WIDTH;
		pb += SSE_WIDTH;
	}
	_mm_storeu_ps( &sum[0], ss );

	for( int i = limit; i < len; i++ )
	{
		sum[0] += a[i] * b[i];
	}

	return sum[0] + sum[1] + sum[2] + sum[3];
}
float Ranf( float low, float high )
{
        float r = (float) rand();               // 0 - RAND_MAX
        float t = r  /  (float) RAND_MAX;       // 0. - 1.

        return   low  +  t * ( high - low );
}

int Ranf( int ilow, int ihigh )
{
        float low = (float)ilow;
        float high = ceil( (float)ihigh );

        return (int) Ranf(low,high);
}
void TimeOfDaySeed( )
{
	unsigned int seed = (unsigned int)( 1000.*3 );    // milliseconds
	srand( seed );
}
int main(){
	TimeOfDaySeed();
	omp_set_num_threads( NUMT );	// same as # of sections
    float *one  = new float [ARRAYSIZE];
	float *two = new float [ARRAYSIZE];
	float *resultArr = new float [ARRAYSIZE];
	float result = 0.;
	float *resultArray = new float [ARRAYSIZE];
	// fill the random-value arrays:
	for( int n = 0; n < ARRAYSIZE; n++ )
	{
		one[n]  = Ranf(  0.0f,  15.f );
		two[n]  = Ranf(  0.0f,  15.f );
	}
double timeResult1 = 9999999999999;
	for (int r= 0; r < 100; r++){
		 double time0 = omp_get_wtime( );
			for( int n = 0; n < ARRAYSIZE; n++ )
			{
				// resultArr[n] = one[n]*two[n];
				// result += resultArr[n];
				result += one[n]*two[n];
				// printf("%5.5f * %5.5f  %5.5f \n",one[n],two[n],result);
			}
			// for( int n = 0; n < ARRAYSIZE; n++ )
			// {
			// 	result += resultArr[n];
			// 	// printf("%5.5f * %5.5f  %5.5f \n",one[n],two[n],result);
			// }
		 double time1 = omp_get_wtime( );
		 if((time1 - time0) < timeResult1){
			 timeResult1 = time1 - time0;
		 }
	}
	// printf("%5.10f,",timeResult1);
	double timeResult2=9999999;
	// 	#pragma omp parallel for
	// for (int r= 0; r < 100; r++){
		float check;
		double time2 = omp_get_wtime( );
		#pragma omp parallel for reduction(+:check)
		for(int i=0; i<omp_get_num_threads(); i=i+ ARRAYSIZE/omp_get_num_threads())
		{
		// #pragma omp parallel
		// float check =SimdMulSum( one+i, two+i,  ARRAYSIZE/omp_get_num_threads() );
		check =SimdMulSum( one+i, two+i,  ARRAYSIZE/omp_get_num_threads() );
		}
		double time3 = omp_get_wtime( );
		 if((time2 - time3) < timeResult2){
			 timeResult2 = time3 - time2;
		 }
	// }
	printf("%d, %d, %5.10f,%5.10f\n",NUMT, ARRAYSIZE,timeResult2,timeResult1);
	// printf("%5.10f %5.10f\n",result, check);
}