/*=========================================================================
 * window functions
 *
 * author: Martin Schiffner
 * date: 2010-10-26
 * modified: 2022-02-18
 *
 *=======================================================================*/

#include "gpu_bf_windows.cuh"

// TODO: add window parameter

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 0.) apply window
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__device__ t_float_gpu apply( const t_window_ptr window, const t_float_gpu l, const t_float_gpu L_over_2 )
{

	// early exit for points outside window
	if( abs( l ) >= L_over_2 ) return 0.0f;

	// return value of window
	return window( l, L_over_2 );

} // __device__ t_float_gpu apply( const t_window_ptr window, const t_float_gpu l, const t_float_gpu L_over_2 )

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 1.) boxcar window
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__device__ t_float_gpu boxcar( const t_float_gpu l, const t_float_gpu L_over_2 )
{

	// equation: 1.0f
	// integral: 2.0f * L_over_2

	// return value of window
// 	return 1.0f / ( 2.0f * L_over_2 );
    return 1.0f;

} // __device__ t_float_gpu boxcar( const t_float_gpu l, const t_float_gpu L_over_2 )

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 2.) Hann (raised-cosine) window
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__device__ t_float_gpu hann( const t_float_gpu l, const t_float_gpu L_over_2 )
{

    // equation: ( 1.0f + cos( M_PI * l / L_over_2 ) ) / 2.0f

    // return value of window
    return ( 1.0f + __cosf( M_PI * l / L_over_2 ) ) / 2.0f;

} // __device__ t_float_gpu hann( const t_float_gpu l, const t_float_gpu L_over_2 )

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 3.) Tukey window
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__device__ t_float_gpu tukey_20( const t_float_gpu l, const t_float_gpu L_over_2 )
{

    // roll-off factor
    t_float_gpu roll_off_factor = 0.2f;

    // auxiliary values
    t_float_gpu positions_over_halfwidth_abs_thresh = 1.0f - roll_off_factor;
    t_float_gpu ratio = abs( l ) / L_over_2;

    // constant part
    if( ratio <= positions_over_halfwidth_abs_thresh ) return 1.0f;
    // assertion: ratio > positions_over_halfwidth_abs_thresh

    // taper part
    return ( 1.0f + __cosf( M_PI * ( ratio - positions_over_halfwidth_abs_thresh ) / roll_off_factor ) ) / 2.0f;

} // __device__ t_float_gpu tukey_20( const t_float_gpu l, const t_float_gpu L_over_2 )

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 4.) cosine window
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__device__ t_float_gpu cosine( const t_float_gpu l, const t_float_gpu L_over_2 )
{

	// equation: cos( M_PI * l / ( 2.0f * L_over_2 ) )
	// integral: 4.0f * L_over_2 / M_PI

	// return value of window
	return __cosf( M_PI * l / ( 2.0f * L_over_2 ) ) / ( 4.0f * L_over_2 / M_PI );

} // __device__ t_float_gpu cosine( const t_float_gpu l, const t_float_gpu L_over_2 )

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 5.) Bartlett (triangular) window
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__device__ t_float_gpu bartlett( const t_float_gpu l, const t_float_gpu L_over_2 )
{

	// equation: 1.0f - abs( l ) / L_over_2
	// integral: L_over_2

	// return value of window
	return ( 1.0f - abs( l ) / L_over_2 ) / L_over_2;

} // __device__ t_float_gpu bartlett( const t_float_gpu l, const t_float_gpu L_over_2 )
