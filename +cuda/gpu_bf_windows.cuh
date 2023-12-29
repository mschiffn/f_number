//-----------------------------------------------------------------------------
// normalized window functions (macros and function prototypes)
//-----------------------------------------------------------------------------
// 0: boxcar
// 1: hann
// 2: tukey_20
// 3: cosine
// 4: bartlett
//
// copy device function pointer array d_windows to host using:
//	t_window_ptr h_windows[ N_WINDOWS ];	// host-side function pointers to __device__ window functions
//	cudaMemcpyFromSymbol( &h_windows, d_windows, N_WINDOWS * sizeof( t_window_ptr ) );
//
#ifndef __GPU_BF_WINDOWS_CUH__
#define __GPU_BF_WINDOWS_CUH__

#include <cuda.h>
#include "gpu_bf_types.cuh"

// number of window functions
#define N_WINDOWS 5

// window configuration
typedef struct
{
	int index;
	double parameter;
} t_window_config;

// window function pointer type
typedef t_float_gpu (*t_window_ptr)( const t_float_gpu, const t_float_gpu, const t_float_gpu );

// device window application function
__device__ t_float_gpu apply( const t_window_ptr window, const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter );

// device window functions
__device__ t_float_gpu boxcar( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter );
__device__ t_float_gpu hann( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter );
__device__ t_float_gpu tukey_20( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter );
__device__ t_float_gpu cosine( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter );
__device__ t_float_gpu bartlett( const t_float_gpu l, const t_float_gpu L_over_2, const t_float_gpu parameter );

// device function pointer array to device window functions
__device__ t_window_ptr d_windows[ N_WINDOWS ] = { boxcar, hann, tukey_20, cosine, bartlett };

#endif
