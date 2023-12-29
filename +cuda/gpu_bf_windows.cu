/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * window functions
 *
 * author: Martin Schiffner
 * date: 2010-10-26
 * modified: 2023-12-28
 *
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "gpu_bf_windows.cuh"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 0.) apply window
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__device__ t_float_gpu apply( const t_window_ptr window, const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )
{

	// early exit for points outside window
	if( abs( l ) >= L_over_2 ) return 0.0f;

	// return value of window
	return window( l, L_over_2, parameter );

} // __device__ t_float_gpu apply( const t_window_ptr window, const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// types of windows
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//----------------------------------------------------------------------------------------------------------------------------------------------
// 0.) boxcar window
//----------------------------------------------------------------------------------------------------------------------------------------------
__device__ t_float_gpu boxcar( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )
{

	// equation: 1.0f
	// integral: 2.0f * L_over_2

	// return value of window
// 	return 1.0f / ( 2.0f * L_over_2 );
    return 1.0f;

} // __device__ t_float_gpu boxcar( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )

//----------------------------------------------------------------------------------------------------------------------------------------------
// 1.) Hann (raised-cosine) window
//----------------------------------------------------------------------------------------------------------------------------------------------
__device__ t_float_gpu hann( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )
{

    // equation: ( 1.0f + cos( M_PI * l / L_over_2 ) ) / 2.0f

    // return value of window
    return ( 1.0f + __cosf( M_PI * l / L_over_2 ) ) / 2.0f;

} // __device__ t_float_gpu hann( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )

//----------------------------------------------------------------------------------------------------------------------------------------------
// 2.) Tukey window
//----------------------------------------------------------------------------------------------------------------------------------------------
__device__ t_float_gpu tukey( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu fraction_cosine )
{

	// auxiliary values
	t_float_gpu positions_over_halfwidth_abs_thresh = 1.0f - fraction_cosine;
	t_float_gpu ratio = abs( l ) / L_over_2;

	// constant part
	if( ratio <= positions_over_halfwidth_abs_thresh ) return 1.0f;
	// assertion: ratio > positions_over_halfwidth_abs_thresh

	// taper part
	return ( 1.0f + __cosf( M_PI * ( ratio - positions_over_halfwidth_abs_thresh ) / fraction_cosine ) ) / 2.0f;

} // __device__ t_float_gpu tukey( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu fraction_cosine )

//----------------------------------------------------------------------------------------------------------------------------------------------
// 3.) cosine window
//----------------------------------------------------------------------------------------------------------------------------------------------
__device__ t_float_gpu cosine( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )
{

	// equation: cos( M_PI * l / ( 2.0f * L_over_2 ) )
	// integral: 4.0f * L_over_2 / M_PI

	// return value of window
	return __cosf( M_PI * l / ( 2.0f * L_over_2 ) ) / ( 4.0f * L_over_2 / M_PI );

} // __device__ t_float_gpu cosine( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )

//----------------------------------------------------------------------------------------------------------------------------------------------
// 4.) Bartlett (triangular) window
//----------------------------------------------------------------------------------------------------------------------------------------------
__device__ t_float_gpu bartlett( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )
{

	// equation: 1.0f - abs( l ) / L_over_2
	// integral: L_over_2

	// return value of window
	return ( 1.0f - abs( l ) / L_over_2 );

} // __device__ t_float_gpu bartlett( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter )
