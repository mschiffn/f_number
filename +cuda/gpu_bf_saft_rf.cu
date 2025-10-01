/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % gpu_bf_saft_rf.cu
 %
 % Form an ultrasound image from
 % a complete synthetic aperture scan by using
 % the delay-and-sum (DAS) algorithm in
 % the Fourier domain.
 %
 % The Fourier domain permits:
 %  1.) independent focusing of different frequencies;
 %  2.) exact corrections of the arrival times independent of the sampling rate;
 %  3.) elimination of out-of-band noise; and
 %  4.) usage of frequency-dependent apodization weights.
 %
 % A uniform linear transducer array is assumed.
 %
 % -----------------------------------------------------------------------------
 % USAGE:
 % -----------------------------------------------------------------------------
 % minimal:
 % image = gpu_bf_saft_rf( positions_x, positions_z, data_RF, f_s, element_width, element_pitch, c_0 );
 %
 % maximal:
 % [ image, weights_tx, weights_rx ] = gpu_bf_saft_rf( positions_x, positions_z, data_RF, f_s, element_width, element_pitch, c_0, f_bounds, index_t0, window_tx, window_rx, F_number_tx, F_number_rx, normalization_tx, normalization_rx, index_gpu, verbosity );
 %
 % -----------------------------------------------------------------------------
 % INPUTS:
 % -----------------------------------------------------------------------------
 % REQUIRED
 %  00.) positions_x:           lateral voxel positions (m)
 %  01.) positions_z:           axial voxel positions (m)
 %  02.) data_RF:               RF data (3d array; 1st dimension: time, 2nd dimension: receive array element, 3rd dimension: transmit array element)
 %  03.) f_s:                   sampling rate of the RF data (Hz)
 %  04.) element_width:         element width of the uniform linear array (m)
 %  05.) element_pitch:         element pitch of the uniform linear array (m)
 %  06.) c_0:                   speed of sound (m/s)
 %
 % OPTIONAL
 %  07.) f_bounds:              frequency bounds (Hz) (1d array; 1st element: lower frequency bound, 2nd element: upper frequency bound)
 %  08.) index_t0:              time index of the sample extracted from the focused RF signal (1)
 %  09.) window_tx:             window function for transmit apodization ( object of class windows.window )
 %  10.) window_rx:             window function for receive apodization ( object of class windows.window )
 %  11.) F_number_tx:           transmit F-number ( object of class f_numbers.f_number )
 %  12.) F_number_rx:           receive F-number ( object of class f_numbers.f_number )
 %  13.) normalization_tx:      normalization of the transmit apodization weights ( object of class normalizations.normalization )
 %  14.) normalization_rx:      normalization of the receive apodization weights ( object of class normalizations.normalization )
 %  15.) index_gpu:             index of the GPU to be used (1)
 %  16.) verbosity:             verbosity level of output (1)
 %
 % -----------------------------------------------------------------------------
 % OUTPUTS:
 % -----------------------------------------------------------------------------
 %  00.) image:                 complex-valued image
 %  01.) weights_tx:            sum of transmit apodization weights for each voxel
 %  02.) weights_rx:            sum of receive apodization weights for each voxel
 %
 % -----------------------------------------------------------------------------
 % ABOUT:
 % -----------------------------------------------------------------------------
 % author: Martin Schiffner
 % date: 2010-12-17
 % modified: 2025-09-09
 % All rights reserved!
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

// CUDA, cuFFT, and cuBLAS
#include <cuda.h>
#include <cufft.h>
#include <math.h>
#include <cublas_v2.h>

#include "gpu_bf_saft_rf.cuh"

// window functions
#include "gpu_bf_windows.cu"

// error handling
#include "gpu_bf_error_handling.cu"

// helper functions
#include "gpu_bf_helper_functions.cu"

#ifdef _CHAR16T
#define CHAR16_T
#endif
#include "mex.h"

#define N_THREADS_X 8
#define N_THREADS_Y 8
#define N_THREADS_PER_BLOCK 64

#define REVISION "1.3"
#define DATE "2025-08-25"
#define BYTES_PER_MEBIBYTE 1048576
#define BYTES_PER_KIBIBYTE 1024

// kernels
#include "das_kernel_distances_sw.cu"			// propagation distances of spherical wave
#include "das_kernel_f_number_constant.cu"		// fixed F-number
#include "saft_kernel_phase_shifts.cu"			// complex-valued apodization weights
#include "das_kernel_normalization.cu"			// normalization of the image

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// MEX gateway function
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 0.) check number of devices supporting CUDA
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// a) check number of devices supporting CUDA
	int N_devices = 0;
	checkCudaRTErrors( cudaGetDeviceCount( &N_devices ) );

	// b) early exit if no CUDA-enabled GPU was detected
	if( N_devices < 1 ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:NoCUDADevices", "Could not find any CUDA-enabled devices!" );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) check and assign required inputs
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// check for correct number of input arguments
	if( nrhs < 7 ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:nrhs", "At least 7 inputs are required." );
	if( nrhs > 17 ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:nrhs", "At most 17 inputs are allowed." );

	// check for correct number of output arguments
	if( nlhs > 3 ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:nlhs", "At most 3 outputs are allowed." );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 0.) lateral voxel positions (m)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with two dimensions
	if( !( mxIsDouble( prhs[ 0 ] ) && !mxIsComplex( prhs[ 0 ] ) && ( mxGetNumberOfDimensions( prhs[ 0 ] ) == 2 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidLateralPositions", "Lateral positions must be a real-valued matrix of type double!" );
	}

	// b) lateral positions
	const double* const positions_x_dbl = ( double * ) mxGetData( prhs[ 0 ] );

	// c) number of positions on x-axis
	const int N_positions_x = mxGetN( prhs[ 0 ] );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 1.) axial voxel positions (m)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with two dimensions
	if( !( mxIsDouble( prhs[ 1 ] ) && !mxIsComplex( prhs[ 1 ] ) && ( mxGetNumberOfDimensions( prhs[ 1 ] ) == 2 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidAxialPositions", "Axial positions must be a real-valued matrix of type double!" );
	}

	// b) axial positions
	const double* const positions_z_dbl = ( double * ) mxGetData( prhs[ 1 ] );

	// c) number of positions on z-axis
	const int N_positions_z = mxGetN( prhs[ 1 ] );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 2.) RF data [ 3d array; 1st dimension: time, 2nd dimension: rx array element index, 3rd dimension: tx array element index ]
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with two dimensions
	if( !( mxIsDouble( prhs[ 2 ] ) && !mxIsComplex( prhs[ 2 ] ) && ( mxGetNumberOfDimensions( prhs[ 2 ] ) == 3 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidRFData", "RF data must be a real-valued array of type double with 3 dimensions!" );
	}

	// b) create pointer to real data in RF data
	const double* const data_RF_dbl = ( double * ) mxGetData( prhs[ 2 ] );

	// c) read dimensions of three-dimensional RF data array
	const mwSize* const RF_data_dim = mxGetDimensions( prhs[ 2 ] );

	const int N_samples_t = RF_data_dim[ 0 ];  // number of samples per RF line in the time domain
	const int N_el_rx = RF_data_dim[ 1 ];      // number of rx array elements
	const int N_el_tx = RF_data_dim[ 2 ];      // number of tx array elements

	// d) ensure equal numbers of transmit and receive
	if( N_el_rx != N_el_tx )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidRFData", "Size of RF data must be equal along 2nd and 3rd dimensions!" );
	}

	const int N_signals = N_el_tx * N_el_rx;  // number of RF signals

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 3.) sampling rate of the RF data (Hz)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 3 ] ) && !mxIsComplex( prhs[ 3 ] ) && mxIsScalar( prhs[ 3 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidSamplingRate", "Sampling rate must be a real-valued scalar of type double!" );
	}

	// b) sampling rate
	const double f_s_dbl = *( ( double * ) mxGetData( prhs[ 3 ] ) );

	// c) ensure positive sampling rate
	if( f_s_dbl <= 0 ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidSamplingRate", "Sampling rate must be positive!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 4.) element width of the uniform linear array (m)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 4 ] ) && !mxIsComplex( prhs[ 4 ] ) && mxIsScalar( prhs[ 4 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidElementWidth", "Element width must be a real-valued scalar of type double!" );
	}

	// b) element width
	const double element_width_dbl = *( ( double * ) mxGetData( prhs[ 4 ] ) );

	// c) ensure validity
	if( !( ( element_width_dbl > 0 ) && mxIsFinite( element_width_dbl ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidElementWidth", "Element width must be positive and finite!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 5.) element pitch of the uniform linear array (m)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 5 ] ) && !mxIsComplex( prhs[ 5 ] ) && mxIsScalar( prhs[ 5 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidElementPitch", "Element pitch must be a real-valued scalar of type double!" );
	}

	// b) element pitch
	const double element_pitch_dbl = *( ( double * ) mxGetData( prhs[ 5 ] ) );

	// c) ensure validity
	if( !( ( element_pitch_dbl > element_width_dbl ) && mxIsFinite( element_pitch_dbl ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidElementPitch", "Element pitch must be finite and larger than element width!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 6.) speed of sound (m/s)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 6 ] ) && !mxIsComplex( prhs[ 6 ] ) && mxIsScalar( prhs[ 6 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidSoundSpeed", "Speed of sound must be a real-valued scalar of type double!" );
	}

	// b) speed of sound
	const double c_0_dbl = *( ( double * ) mxGetData( prhs[ 6 ] ) );

	// c) ensure validity
	if( !( ( c_0_dbl > 0 ) && mxIsFinite( c_0_dbl ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidSoundSpeed", "Average speed of sound must be positive and finite!" );
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 2.) initialize and assign optional inputs
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 7.) frequency bounds (Hz) (1d array; 1st element: lower frequency bound, 2nd element: upper frequency bound)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default values
	double f_lb_dbl = DBL_MIN;		// lower bound: smallest double greater than zero
	double f_ub_dbl = f_s_dbl / 2; 	// upper bound: Nyquist frequency

	// b) check presence of nonempty argument
	if( nrhs >= 8 && !mxIsEmpty( prhs[ 7 ] ) )
	{
		// ensure real-valued mxDOUBLE_CLASS with two elements
		if( !( mxIsDouble( prhs[ 7 ] ) && !mxIsComplex( prhs[ 7 ] ) && ( mxGetNumberOfElements( prhs[ 7 ] ) == 2 ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyBounds", "Frequency bounds must be a real-valued vector with two entries of type double!" );
		}

		// frequency bounds
		const double* const f_bounds_dbl = ( double * ) mxGetData( prhs[ 7 ] );

		// extract frequency bounds
		f_lb_dbl = f_bounds_dbl[ 0 ];
		f_ub_dbl = f_bounds_dbl[ 1 ];
	}

	// c) ensure positive lower bound and valid upper bound
	if( f_lb_dbl <= 0 ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyBound", "Lower frequency bound must be positive!" );
	if( f_ub_dbl <= f_lb_dbl ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyBound", "Upper frequency bound must exceed the lower frequency bound!" );
	if( f_ub_dbl > f_s_dbl / 2 ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyBound", "Upper frequency bound must not exceed the Nyquist frequency!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 8.) time index of the sample extracted from the focused RF signal
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	double index_t0_dbl = 0.0;

	// b) check specified value
	if( nrhs >= 9 && !mxIsEmpty( prhs[ 8 ] ) )
	{
		// ensure real-valued mxDOUBLE_CLASS with one element
		if( !( mxIsDouble( prhs[ 8 ] ) && !mxIsComplex( prhs[ 8 ] ) && mxIsScalar( prhs[ 8 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidTimeShift", "Time index must be a real-valued scalar of type double!" );
		}

		// read value
		index_t0_dbl = *( ( double * ) mxGetData( prhs[ 8 ] ) );
	}

	// c) ensure nonnegative and finite time index
	if( !( ( index_t0_dbl >= 0 ) && mxIsFinite( index_t0_dbl ) ) ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidTimeShift", "Time index must be nonnegative!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 9.) window function for the transmit apodization ( object of class windows.window )
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	t_window_config window_tx_config = { .index = 0, .parameter = 0.0 };

	// b) check specified value
	if( nrhs >= 10 && !mxIsEmpty( prhs[ 9 ] ) )
	{
		// ensure single object of class windows.window
		if( !( isa( prhs[ 9 ], "windows.window" ) && mxIsScalar( prhs[ 9 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Transmit window must be a single object of the class windows.window!");
		}

		// map classes to window id
		if( get_window_id( &window_tx_config, prhs[ 9 ] ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Could not read parameters of transmit window!" );
		}
	}

	// c) check validity of window configuration
	if( !( ( window_tx_config.index >= 0 ) && ( window_tx_config.index < N_WINDOWS ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Unknown transmit window function!" );
	}
	if( !( ( window_tx_config.parameter >= 0 ) && ( window_tx_config.parameter <= 1 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Invalid transmit window parameter!" );
	}

	// d) name of window function
	char window_tx_name[ 128 ];
	if( string( prhs[ 9 ], window_tx_name, 128 * sizeof( char ) ) ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Error reading transmit window name!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 10.) window function for the receive apodization ( object of class windows.window )
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	t_window_config window_rx_config = { .index = 0, .parameter = 0.0 };

	// b) check specified value
	if( nrhs >= 11 && !mxIsEmpty( prhs[ 10 ] ) )
	{
		// ensure single object of class windows.window
		if( !( isa( prhs[ 10 ], "windows.window" ) && mxIsScalar( prhs[ 10 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Receive window must be a single object of the class windows.window!");
		}

		// map classes to window id
		if( get_window_id( &window_rx_config, prhs[ 10 ] ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Could not read parameters of receive window!" );
		}
	}

	// c) check validity of window configuration
	if( !( ( window_rx_config.index >= 0 ) && ( window_rx_config.index < N_WINDOWS ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Unknown receive window function!" );
	}
	if( !( ( window_rx_config.parameter >= 0 ) && ( window_rx_config.parameter <= 1 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Invalid receive window parameter!" );
	}

	// d) name of window function
	char window_rx_name[ 128 ];
	if( string( prhs[ 10 ], window_rx_name, 128 * sizeof( char ) ) ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Error reading receive window name!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 11.) transmit F-number ( object of class f_numbers.f_number )
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) create default transmit F-number
	mxArray* angles_lb_deg_matlab = mxCreateDoubleScalar( 70.0 );
	mxArray* F_numbers_ub_matlab = mxCreateDoubleScalar( 1.5 );
	mxArray* distances_deg_matlab = mxCreateDoubleScalar( 10.0 );
	mxArray* rhs[ 3 ] = { angles_lb_deg_matlab, F_numbers_ub_matlab, distances_deg_matlab };
	mxArray* f_number_tx_matlab = NULL;
	mexCallMATLAB( 1, &f_number_tx_matlab, 3, rhs, "f_numbers.grating.angle_lb" );
	bool f_number_tx_constant = false;

	// b) check specified value
	if( nrhs >= 12 && !mxIsEmpty( prhs[ 11 ] ) )
	{
		// ensure single object of class f_numbers.f_number
		if( !( isa( prhs[ 11 ], "f_numbers.f_number" ) && mxIsScalar( prhs[ 11 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFNumber", "Transmit F-number must be a single object of the class f_numbers.f_number!");
		}

		// destroy default transmit F-number
		// mxDestroyArray( f_number_tx_matlab );

		// make deep copy of transmit F-number to avoid problems w/ const modifier
		f_number_tx_matlab = mxDuplicateArray( prhs[ 11 ] );
	}

	// c) check if specified transmit F-number is fixed
	if( isa( f_number_tx_matlab, "f_numbers.constant" ) )
	{
		f_number_tx_constant = true;
	}

	// name of transmit F-number
	char f_number_tx_name[ 128 ];
	if( string( f_number_tx_matlab, f_number_tx_name, 128 * sizeof(char) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFNumber", "Error reading transmit F-number name!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 12.) receive F-number ( object of class f_numbers.f_number )
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) create default receive F-number
	mxArray* f_number_rx_matlab = NULL;
	mexCallMATLAB( 1, &f_number_rx_matlab, 3, rhs, "f_numbers.grating.angle_lb" );
	bool f_number_rx_constant = false;

	// b) check specified value
	if( nrhs >= 13 && !mxIsEmpty( prhs[ 12 ] ) )
	{
		// ensure single object of class f_numbers.f_number
		if( !( isa( prhs[ 12 ], "f_numbers.f_number" ) && mxIsScalar( prhs[ 12 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFNumber", "Receive F-number must be a single object of the class f_numbers.f_number!");
		}

		// destroy default receive F-number
		// mxDestroyArray( f_number_rx_matlab );

		// make deep copy of receive F-number to avoid problems w/ const modifier
		f_number_rx_matlab = mxDuplicateArray( prhs[ 12 ] );
	}

	// c) check if specified receive F-number is fixed
	if( isa( f_number_rx_matlab, "f_numbers.constant" ) )
	{
		f_number_rx_constant = true;
	}

	// name of receive F-number
	char f_number_rx_name[ 128 ];
	if( string( f_number_rx_matlab, f_number_rx_name, 128 * sizeof(char) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFNumber", "Error reading receive F-number name!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 13.) normalization of the transmit apodization weights ( object of class normalizations.normalization )
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default values
	bool normalization_tx = true;
	t_window_config window_f_tx_config = { .index = 0, .parameter = 0.0 };
	char normalization_tx_name[ 128 ] = { "on" };

	// b) check specified value
	if( nrhs >= 14 && !mxIsEmpty( prhs[ 13 ] ) )
	{
		// ensure scalar object of class normalizations.normalization
		if( !( isa( prhs[ 13 ], "normalizations.normalization" ) && mxIsScalar( prhs[ 13 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidNormalization", "Transmit normalization must be a single object of the class normalizations.normalization!");
		}

		// name of transmit normalization
		if( string( prhs[ 13 ], normalization_tx_name, 128 * sizeof(char) ) ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidNormalization", "Error reading transmit normalization name!" );

		// deactivate transmit normalization
		if( mxIsClass( prhs[ 13 ], "normalizations.off" ) ) normalization_tx = false;

		// map classes to window id
		if( mxIsClass( prhs[ 13 ], "normalizations.windowed" ) )
		{
			if( get_window_id( &window_f_tx_config, mxGetProperty( prhs[ 13 ], 0, "window" ) ) )
			{
				mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyWindow", "Could not read parameters of window for transmit normalization!" );
			}
		}
	}

	// c) check validity of window configuration
	if( !( ( window_f_tx_config.index >= 0 ) && ( window_f_tx_config.index < N_WINDOWS ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyWindow", "Unknown frequency window function!" );
	}
	if( !( ( window_f_tx_config.parameter >= 0 ) && ( window_f_tx_config.parameter <= 1 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Invalid window parameter!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 14.) normalization of the receive apodization weights ( object of class normalizations.normalization )
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default values
	bool normalization_rx = true;
	t_window_config window_f_rx_config = { .index = 0, .parameter = 0.0 };
	char normalization_rx_name[ 128 ] = { "on" };

	// b) check specified value
	if( nrhs >= 15 && !mxIsEmpty( prhs[ 14 ] ) )
	{
		// ensure scalar object of class normalizations.normalization
		if( !( isa( prhs[ 14 ], "normalizations.normalization" ) && mxIsScalar( prhs[ 14 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidNormalization", "Receive normalization must be a single object of the class normalizations.normalization!");
		}

		// name of receive normalization
		if( string( prhs[ 14 ], normalization_rx_name, 128 * sizeof(char) ) ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidNormalization", "Error reading receive normalization name!" );

		// deactivate receive normalization
		if( mxIsClass( prhs[ 14 ], "normalizations.off" ) ) normalization_rx = false;

		// map classes to window id
		if( mxIsClass( prhs[ 14 ], "normalizations.windowed" ) )
		{
			if( get_window_id( &window_f_rx_config, mxGetProperty( prhs[ 14 ], 0, "window" ) ) )
			{
				mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyWindow", "Could not read parameters of window for receive normalization!" );
			}
		}
	}

	// c) check validity of window configuration
	if( !( ( window_f_rx_config.index >= 0 ) && ( window_f_rx_config.index < N_WINDOWS ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyWindow", "Unknown frequency window function!" );
	}
	if( !( ( window_f_rx_config.parameter >= 0 ) && ( window_f_rx_config.parameter <= 1 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidWindow", "Invalid window parameter!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 15.) index of the GPU
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	int index_gpu = 0;

	// b) check specified value
	if( nrhs >= 16 && !mxIsEmpty( prhs[ 15 ] ) )
	{
		// ensure real-valued mxDOUBLE_CLASS with one element
		if( !( mxIsDouble( prhs[ 15 ] ) && !mxIsComplex( prhs[ 15 ] ) && mxIsScalar( prhs[ 15 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidCUDADevice", "Device index must be a real-valued scalar of type double!" );
		}

		// read device index
		index_gpu = (int) *( ( double * ) mxGetData( prhs[ 15 ] ) );
	}

	// c) early exit if selected GPU does not exist
	if( !( ( index_gpu >= 0 ) || ( index_gpu < N_devices ) ) ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidCUDADevice", "Selected device does not exist!" );

	// d) select device
	checkCudaRTErrors( cudaSetDevice( index_gpu ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 16.) verbosity level of output
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	char verbosity = 1;

	// b) check specified value
	if( nrhs >= 17 && !mxIsEmpty( prhs[ 16 ] ) )
	{
		// ensure real-valued mxDOUBLE_CLASS with one element
		if( !( mxIsDouble( prhs[ 16 ] ) && !mxIsComplex( prhs[ 16 ] ) && mxIsScalar( prhs[ 16 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidVerbosityLevel", "Verbosity level must be a real-valued scalar of type double!" );
		}

		// read verbosity level
		verbosity = (char) *( ( double * ) mxGetData( prhs[ 16 ] ) );
	}

	// c) early exit if selected verbosity level does not exist
	if( ( verbosity < 0 ) || ( verbosity >= 4 ) ) mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidVerbosityLevel", "Selected verbosity level does not exist!" );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 3.) compute dependent parameters and display status information
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) array geometry
	//------------------------------------------------------------------------------------------------------------------------------------------
	// half-number of receive elements
	const double M_el_rx = ( N_el_rx - 1.0 ) / 2.0;

	// compute lateral positions of array elements
	float* pos_rx_ctr_x = ( float * ) mxMalloc( sizeof( float ) * N_el_rx );
	float* pos_rx_lb_x = ( float * ) mxMalloc( sizeof( float ) * N_el_rx );
	float* pos_rx_ub_x = ( float * ) mxMalloc( sizeof( float ) * N_el_rx );

	for( int index_element = 0; index_element < N_el_rx; index_element++ )
	{
		pos_rx_ctr_x[ index_element ] = (float) ( element_pitch_dbl * ( index_element - M_el_rx ) );
		pos_rx_lb_x[ index_element ] = pos_rx_ctr_x[ index_element ] - element_width_dbl / 2;
		pos_rx_ub_x[ index_element ] = pos_rx_ctr_x[ index_element ] + element_width_dbl / 2;
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// b) compute order of DFT; employ zero-padding to achieve odd order
	//------------------------------------------------------------------------------------------------------------------------------------------
	// maximum relative time shift in the electronic focusing ( zero padding in DFT )
// TODO: correct for SAFT?
	const int N_samples_t_pad = ceil( ( sqrt( ( pos_rx_ctr_x[ N_el_rx - 1 ] - pos_rx_ctr_x[ 0 ] ) * ( pos_rx_ctr_x[ N_el_rx - 1 ] - pos_rx_ctr_x[ 0 ] ) + positions_z_dbl[ 1 ] * positions_z_dbl[ 1 ] ) - positions_z_dbl[ 1 ] ) * f_s_dbl / c_0_dbl );

	// number of points in DFT
	const int N_order_dft = N_samples_t + N_samples_t_pad + 1 - ( ( N_samples_t + N_samples_t_pad ) % 2 );
	// assertion: N_order_dft is odd

	if( verbosity > 1 ) mexPrintf( "N_order_dft = %d\n", N_order_dft );

	// number of distinct samples in DFT of real-valued signal
	const int N_samples_dft = N_order_dft / 2 + 1;

	// arguments of complex exponential for focusing
	const double argument_coefficient_dbl = 2 * M_PI * f_s_dbl / ( N_order_dft * c_0_dbl );
	const double argument_add_dbl = 2 * M_PI * index_t0_dbl / N_order_dft;

	//------------------------------------------------------------------------------------------------------------------------------------------
	// c) compute parameters of bandpass filter
	//------------------------------------------------------------------------------------------------------------------------------------------
	// lower frequency index in the RF signals
	const int index_f_lb = ceil( f_lb_dbl * N_order_dft / f_s_dbl );

	// upper frequency index in the RF signals
	const int index_f_ub = floor( f_ub_dbl * N_order_dft / f_s_dbl );

	// check validity of frequency indices
	if( !( ( index_f_lb > 0 ) && ( index_f_lb <= index_f_ub ) && ( index_f_ub < N_samples_dft ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidFrequencyBounds", "Frequency bounds lead to invalid indices!" );
	}

	// number of non-redundant samples of DFT of order N_order_dft after bandpass filter
	const int N_samples_dft_bp  = index_f_ub - index_f_lb + 1;

	//------------------------------------------------------------------------------------------------------------------------------------------
	// d) compute frequency-dependent F-numbers [ in MATLAB ]
	//------------------------------------------------------------------------------------------------------------------------------------------
	// create normalized element pitches in MATLAB workspace
	mxArray* element_pitch_norm_matlab = mxCreateDoubleMatrix( N_samples_dft_bp, 1, mxREAL );

	// frequency axis
	double* element_pitch_norm_dbl = (double *) mxGetData( element_pitch_norm_matlab );

	// iterate frequencies
	for( int index_f = 0; index_f < N_samples_dft_bp; index_f++ )
	{
		// normalized element pitch
		element_pitch_norm_dbl[ index_f ] = element_pitch_dbl * ( index_f_lb + index_f ) * f_s_dbl / ( N_order_dft * c_0_dbl );
	}

	// compute values of the frequency-dependent F-number
	const double* const f_number_tx_values_dbl = compute_values( f_number_tx_matlab, element_pitch_norm_matlab );
	const double* const f_number_rx_values_dbl = compute_values( f_number_rx_matlab, element_pitch_norm_matlab );

	// ensure valid F-numbers
	for( int index_f = 0; index_f < N_samples_dft_bp; index_f++ )
	{
		// ensure valid transmit F-number
		if( !( ( f_number_tx_values_dbl[ index_f ] >= 0 ) && ( f_number_tx_values_dbl[ index_f ] <= DBL_MAX ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidTXFNumber", "Values of transmit F-number must be nonnegative!" );
		}
		// ensure valid receive F-number
		if( !( ( f_number_rx_values_dbl[ index_f ] >= 0 ) && ( f_number_rx_values_dbl[ index_f ] <= DBL_MAX ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft_rf:InvalidRXFNumber", "Values of transmit F-number must be nonnegative!" );
		}
		// mexPrintf( "f_number_tx_values_dbl[%d] = %.2f\n", index_f, f_number_tx_values_dbl[ index_f ] );
		// mexPrintf( "f_number_rx_values_dbl[%d] = %.2f\n", index_f, f_number_rx_values_dbl[ index_f ] );
	}

	// destroy normalized element pitches in MATLAB workspace
	// mxDestroyArray( element_pitch_norm_matlab );
	// element_pitch_norm_dbl = NULL;

	// destroy F-number
	// mxDestroyArray( f_number_rx_matlab );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// e) compute parallelization settings
	//------------------------------------------------------------------------------------------------------------------------------------------
	// number of blocks to process in parallel
	const int N_blocks_x = ceil( ( (double) N_positions_x ) / N_THREADS_X );
	const int N_blocks_z = ceil( ( (double) N_positions_z ) / N_THREADS_Y );
	const int N_blocks_rx = ceil( ( (double) N_el_rx ) / N_THREADS_PER_BLOCK );

	// determine relevant range of frequency indices
	const int N_blocks_Omega_bp = ceil( ((double) N_samples_dft_bp) / N_THREADS_PER_BLOCK );  // number of frequency blocks

	//------------------------------------------------------------------------------------------------------------------------------------------
	// f) print global status information
	//------------------------------------------------------------------------------------------------------------------------------------------
	if( verbosity > 0 )
	{
		mexPrintf( " %s\n", "================================================================================" );
		mexPrintf( " GPU BF Toolbox: DAS, SAFT, v.%s (%s), Status Information\n", REVISION, DATE );
		mexPrintf( " (c) 2010-2025, M. F. Schiffner. All rights reserved.\n" );
		mexPrintf( " %s\n", "================================================================================" );
		mexPrintf( "  %-20s: % 5d x %-5d %5s", "lattice geometry", N_positions_x, N_positions_z, "" );
		mexPrintf( " %-22s: %-8d\n", "number of samples", N_samples_t );
		mexPrintf( "  %-20s: %-11d %7s", "number transmits", N_el_tx, "" );
		mexPrintf( " %-22s: %5.1f MHz\n", "sampling rate", f_s_dbl / 1000000 );
		mexPrintf( "  %-20s: %-11d %7s", "number rx channels", N_el_rx, "" );
		mexPrintf( " %-22s: %5.1f MHz (%d)\n", "lower frequency bound", f_lb_dbl / 1000000, index_f_lb );
		mexPrintf( "  %-20s: %7.2f m/s %7s", "average sound speed", c_0_dbl, "" );
		mexPrintf( " %-22s: %5.1f MHz (%d)\n", "upper frequency bound", f_ub_dbl / 1000000, index_f_ub );
		mexPrintf( "  %-20s: %-6.2f %12s", "additional shift", index_t0_dbl, "" );
		mexPrintf( " %-22s: %-6d\n", "number of frequencies", N_samples_dft_bp );
		mexPrintf( "  %-20s: %s (%d; %3.2f)\n", "window (tx)", window_tx_name, window_tx_config.index, window_tx_config.parameter );
		mexPrintf( "  %-20s: %s (%d; %3.2f)\n", "window (rx)", window_rx_name, window_rx_config.index, window_rx_config.parameter );
		mexPrintf( "  %-20s: %s\n", "F-number (tx)", f_number_tx_name );
		mexPrintf( "  %-20s: %s\n", "F-number (rx)", f_number_rx_name );
		mexPrintf( "  %-20s: %s", "normalization (tx)", normalization_tx_name );
		mexPrintf( " %-22s: %s\n", "normalization (rx)", normalization_rx_name );

		//mexPrintf(" %-22s: %-6.2f um\n", "maximum delta_x", delta_x_max * 1e6);
		//mexPrintf(" %-22s: %-6.2f um\n", "maximum delta_z", delta_z_max * 1e6);
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 4.) infer and display properties of selected device
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) hardware properties
	//------------------------------------------------------------------------------------------------------------------------------------------
	// check properties of chosen GPU
	cudaDeviceProp deviceProp;
	checkCudaRTErrors( cudaGetDeviceProperties( &deviceProp, index_gpu ) );

	// check capabilities for double precision
	// int capability_double = 0;
	// if( ( deviceProp.major == 1 ) && ( deviceProp.minor >= 3 ) ) capability_double = 1;
	// if( deviceProp.major >= 2 ) capability_double = 1;

	//------------------------------------------------------------------------------------------------------------------------------------------
	// b) software properties
	//------------------------------------------------------------------------------------------------------------------------------------------
	// get driver version
	int driverVersion = 0;
	checkCudaRTErrors( cudaDriverGetVersion( &driverVersion ) );

	// get runtime version
	int runtimeVersion = 0;
	checkCudaRTErrors( cudaRuntimeGetVersion( &runtimeVersion ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// c) print GPU status information
	//------------------------------------------------------------------------------------------------------------------------------------------
	if( verbosity > 0 )
	{
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf(" Selected device %-1d of %-1d: %-56s\n", index_gpu, N_devices, deviceProp.name );
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		// mexPrintf("  %-20s: %-19s", "name of device", deviceProp.name, "");
		mexPrintf("  %-20s: %-1d.%-1d %15s", "compute capability", deviceProp.major, deviceProp.minor, "");
		mexPrintf(" %-22s: %-3d\n", "num. multiprocessors", deviceProp.multiProcessorCount);
		// mexPrintf(" %-22s: %-8d\n", "double precision", capability_double);
		mexPrintf("  %-20s: %-8.2f MiByte %3s", "total global memory", ((double) deviceProp.totalGlobalMem) / BYTES_PER_MEBIBYTE, "");
		mexPrintf(" %-22s: %-6.2f KiByte\n", "total constant memory", ((double) deviceProp.totalConstMem) / BYTES_PER_KIBIBYTE);
		mexPrintf("  %-20s: %-6.2f KiByte %5s", "shared mem. / block", ((double) deviceProp.sharedMemPerBlock) / BYTES_PER_KIBIBYTE, "");
		mexPrintf(" %-22s: %-8d\n", "registers per block", deviceProp.regsPerBlock);
		mexPrintf("  %-20s: %2d.%1d %14s", "driver version", driverVersion / 1000, (driverVersion % 100) / 10, "");
		mexPrintf(" %-22s: %2d.%1d\n", "runtime version", runtimeVersion / 1000, (runtimeVersion % 100) / 10);
		mexPrintf("  %-20s: %-8s %10s", "selected precision", "single", "");
		mexPrintf(" %-22s: %-8d\n", "warp size", deviceProp.warpSize);

		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf(" parallelization parameters\n");
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf("  %-20s: %-4d x %-4d %7s", "number of threads", N_THREADS_X, N_THREADS_Y, "");
		mexPrintf(" %-22s: %-4d (%-4d)\n", "threads per blocks", N_THREADS_PER_BLOCK, deviceProp.maxThreadsPerBlock);
		mexPrintf("  %-20s: %-4d x %-4d %7s", "number of blocks", N_blocks_x, N_blocks_z, "");
		mexPrintf(" %-22s: %-4d\n", "frequency blocks", N_blocks_Omega_bp);
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 5.) allocate and initialize device memory / convert data to float and copy to the device
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// size of conversion buffer: use maximum of
		// a) d_positions_x:
		//    sizeof( t_float_gpu ) * N_positions_x
		// b) d_positions_z:
		//    sizeof( t_float_gpu ) * N_positions_z;
		// c) d_data_RF:
		//    sizeof( cufftComplex ) * N_samples_dft * N_el_rx
		// d) d_pos_rx_ctr_x, d_pos_rx_lb_x, d_pos_rx_ub_x:
		//    sizeof( t_float_gpu ) * N_el_rx
		// f) d_image:
		//    sizeof( cufftComplex ) * N_positions_x * N_positions_z
		// g) d_weights_rx:
		//    sizeof( t_float_gpu ) * N_positions_x * N_positions_z
		//
		// => c) > d)
		// => f) > g)
		// => a), d) vs. c)
		// N_f_unique * max( N_points, N_objects ) vs. N_points_occupied * N_objects

		// b) d_indices_grid_FOV_shift[ index_element ]: N_points_occupied * sizeof( int ) [ indices of shifted grid points for each array element ]
// reinterpret the memory
	// const size_t size_bytes_buffer = sizeof( t_float_gpu ) * N_samples_t * N_el_rx;
	// if( DEBUG_MODE ) mexPrintf( "size_bytes_buffer = %.2f MiB (%zu B)\n", ( ( double ) size_bytes_buffer ) / BYTES_PER_MEBIBYTE, size_bytes_buffer );

	// allocate page-locked memory
	// t_float_gpu* buffer = NULL;
	// checkCudaErrors( cudaHostAlloc( (void**) &buffer, size_bytes_buffer, cudaHostAllocDefault ) );

	// create events for time-measurement
	cudaEvent_t start, stop;
	checkCudaRTErrors( cudaEventCreate( &start ) );
	checkCudaRTErrors( cudaEventCreate( &stop ) );

	// place start event into the default stream
	checkCudaRTErrors( cudaEventRecord( start, 0 ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 0.) lateral voxel positions
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes (FLOAT)
	const size_t size_bytes_positions_x = sizeof( t_float_gpu ) * N_positions_x;

	// b) allocate memory
	t_float_gpu* positions_x = (t_float_gpu *) mxMalloc( size_bytes_positions_x );

	// c) convert lateral image coordinates
	for( int l_x = 0; l_x < N_positions_x; l_x++ )
	{
		positions_x[ l_x ] = (t_float_gpu) positions_x_dbl[ l_x ];
	}

	// allocate device memory
	t_float_gpu* d_positions_x = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_positions_x, size_bytes_positions_x ) );
	// checkCudaErrors( cudaMallocPitch( (void **) &d_h_ref_float_complex, &pitch_h_ref, N_f_unique * sizeof( t_float_complex_gpu ), N_points ) );

	// copy data
	checkCudaRTErrors( cudaMemcpy( d_positions_x, positions_x, size_bytes_positions_x, cudaMemcpyHostToDevice ) );
	// checkCudaErrors( cudaMemcpy2D( d_h_ref_float_complex, pitch_h_ref, buffer, N_f_unique * sizeof( t_float_complex_gpu ), N_f_unique * sizeof( t_float_complex_gpu ), N_points, cudaMemcpyHostToDevice ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 1.) axial voxel positions
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes (FLOAT)
	const size_t size_bytes_positions_z = sizeof( t_float_gpu ) * N_positions_z;

	// b) allocate memory
	t_float_gpu* positions_z = (t_float_gpu *) mxMalloc( size_bytes_positions_z );

	// c) conversion to selected data type
	for( int l_z = 0; l_z < N_positions_z; l_z++ )
	{
		positions_z[ l_z ] = (t_float_gpu) positions_z_dbl[ l_z ];
	}

	// d) allocate device memory
	t_float_gpu* d_positions_z = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_positions_z, size_bytes_positions_z ) );

	// e) copy data
	checkCudaRTErrors( cudaMemcpy( d_positions_z, positions_z, size_bytes_positions_z, cudaMemcpyHostToDevice ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 2.) RF data
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes (FLOAT)
	const size_t size_bytes_RF_data = sizeof( cufftComplex ) * N_samples_dft * N_signals;
	const size_t size_bytes_RF_data_host = sizeof( t_float_gpu ) * N_samples_t * N_signals;

	// b) allocate memory
	t_float_gpu* data_RF = (t_float_gpu *) mxMalloc( size_bytes_RF_data_host );

	// c) conversion to selected data type
	for( int k_tx = 0; k_tx < N_el_tx; k_tx++ )
	{
		for( int k_rx = 0; k_rx < N_el_rx; k_rx++ )
		{
			for( int l_z = 0; l_z < N_samples_t; l_z++ )
			{
				data_RF[ ( k_tx * N_el_rx + k_rx ) * N_samples_t + l_z ] = (t_float_gpu) data_RF_dbl[ ( k_tx * N_el_rx + k_rx ) * N_samples_t + l_z ];
			}
		}
	}

	// d) allocate device memory for RF data and its DFT [additional memory for in-place DFT]
	t_float_gpu* d_data_RF = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_data_RF, size_bytes_RF_data ) );

	// e) initialize RF data (necessary because of zero-padding)
	checkCudaRTErrors( cudaMemset( d_data_RF, 0, size_bytes_RF_data ) );

	// f) copy data [correctly pad RF data to allow in place transform]
	const size_t size_bytes_signal = sizeof( t_float_gpu ) * N_samples_t;
	for( int k_tx = 0; k_tx < N_el_tx; k_tx++ )
	{
		for( int k_rx = 0; k_rx < N_el_rx; k_rx++ )
		{
			checkCudaRTErrors( cudaMemcpy( d_data_RF + ( k_tx * N_el_rx + k_rx ) * N_samples_dft * 2, &( data_RF[ ( k_tx * N_el_rx + k_rx ) * N_samples_t ] ), size_bytes_signal, cudaMemcpyHostToDevice ) );
		}
	}
	// N_order_dft / 2 + 1 | data_RF + k_rx * N_samples_t

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 3.) array geometry
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes
	const size_t size_bytes_positions_rx = sizeof( t_float_gpu ) * N_el_rx;

	// b) allocate device memory for array geometry
	t_float_gpu* d_pos_rx_ctr_x = NULL; t_float_gpu* d_pos_rx_lb_x = NULL; t_float_gpu* d_pos_rx_ub_x = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_ctr_x, size_bytes_positions_rx) );
	checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_lb_x, size_bytes_positions_rx) );
	checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_ub_x, size_bytes_positions_rx) );

	// c) copy data
	checkCudaRTErrors( cudaMemcpy( d_pos_rx_ctr_x, pos_rx_ctr_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );
	checkCudaRTErrors( cudaMemcpy( d_pos_rx_lb_x, pos_rx_lb_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );
	checkCudaRTErrors( cudaMemcpy( d_pos_rx_ub_x, pos_rx_ub_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 4.) propagation distances of incident wave
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes
	const size_t size_bytes_distances_tx = sizeof( t_float_gpu ) * N_positions_x * N_positions_z;

	// b) allocate device memory
	t_float_gpu* d_distances_tx = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_distances_tx, size_bytes_distances_tx ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 5.) image
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes
	const size_t size_bytes_DAS_image_RF_complex = sizeof( cufftComplex ) * N_positions_x * N_positions_z;

	// b) allocate device memory
	cufftComplex* d_image = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_image, size_bytes_DAS_image_RF_complex ) );

	// c) initialize complex-valued DAS image
	checkCudaRTErrors( cudaMemset( d_image, 0, size_bytes_DAS_image_RF_complex ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 6.) weights
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes
	size_t size_bytes_weights = sizeof( t_float_gpu ) * N_positions_x * N_positions_z;

	// b) allocate device memory
	t_float_gpu* d_weights_tx = NULL; t_float_gpu* d_weights_rx = NULL;
	checkCudaRTErrors( cudaMalloc( (void**) &d_weights_tx, size_bytes_weights ) );
	checkCudaRTErrors( cudaMalloc( (void**) &d_weights_rx, size_bytes_weights ) );

	// c) initialize weights w/ zeros
	checkCudaRTErrors( cudaMemset( d_weights_tx, 0, size_bytes_weights ) );
	checkCudaRTErrors( cudaMemset( d_weights_rx, 0, size_bytes_weights ) );

	// place stop event into the default stream
	checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
	checkCudaRTErrors( cudaEventSynchronize( stop ) );

	// compute elapsed time for memory transfer host to device (ms)
	float time_transfer_to_device = -1.0f;
	checkCudaRTErrors( cudaEventElapsedTime( &time_transfer_to_device, start, stop ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 7.) memory status
	//------------------------------------------------------------------------------------------------------------------------------------------
	// array geometry
	const size_t size_bytes_geometry_array = 3 * size_bytes_positions_rx;

	// matrix of phase shifts
	const size_t size_bytes_phase = sizeof( cufftComplex ) * N_positions_x * N_positions_z * N_el_rx;

	// total memory consumption
	const size_t size_bytes_total = size_bytes_positions_x + size_bytes_positions_z + size_bytes_RF_data + size_bytes_geometry_array + size_bytes_distances_tx + size_bytes_phase + size_bytes_DAS_image_RF_complex + size_bytes_weights;

	// get available amount of global memory on GPU device
	size_t mem_available = 0;
	checkCudaRTErrors( cudaMemGetInfo( &mem_available, &(deviceProp.totalGlobalMem) ) );

	if( verbosity > 0 )
	{
		mexPrintf( " %s\n", "--------------------------------------------------------------------------------" );
		mexPrintf( " memory consumption (available global memory: %-6.2f MiByte (%6.2f %%))\n", ((double) mem_available) / BYTES_PER_MEBIBYTE, ((double) mem_available) * 100 / deviceProp.totalGlobalMem );
		mexPrintf( " %s\n", "--------------------------------------------------------------------------------" );
		mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "lateral voxel positions", ((double) size_bytes_positions_x) / BYTES_PER_KIBIBYTE, ((double) size_bytes_positions_x) * 100 / size_bytes_total );
		mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "axial voxel positions", ((double) size_bytes_positions_z) / BYTES_PER_KIBIBYTE, ((double) size_bytes_positions_z) * 100 / size_bytes_total );
		mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "RF signals", ((double) size_bytes_RF_data) / BYTES_PER_KIBIBYTE, ((double) size_bytes_RF_data) * 100 / size_bytes_total );
		mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "array geometry", ((double) size_bytes_geometry_array) / BYTES_PER_KIBIBYTE, ((double) size_bytes_geometry_array) * 100 / size_bytes_total );
		mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "TX distances", ((double) size_bytes_distances_tx) / BYTES_PER_KIBIBYTE, ((double) size_bytes_distances_tx) * 100 / size_bytes_total );
		if( !f_number_rx_constant ) mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "phase matrix", ( (double) size_bytes_phase ) / BYTES_PER_KIBIBYTE, ((double) size_bytes_phase) * 100 / size_bytes_total );
		if( !f_number_rx_constant ) mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "DAS image (monofrequent)", ((double) size_bytes_DAS_image_RF_complex) / BYTES_PER_KIBIBYTE, ((double) size_bytes_DAS_image_RF_complex) * 100 / size_bytes_total );
		mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "DAS image", ((double) size_bytes_DAS_image_RF_complex) / BYTES_PER_KIBIBYTE, ((double) size_bytes_DAS_image_RF_complex) * 100 / size_bytes_total );
		mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "weights", ((double) size_bytes_weights) / BYTES_PER_KIBIBYTE, ((double) size_bytes_weights) * 100 / size_bytes_total );
		mexPrintf( "%4s %s  %s\n", "", "-----------------------------", "----------------------------" );
		mexPrintf( "%4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "total", ((double) size_bytes_total) / BYTES_PER_KIBIBYTE, 100.0 );
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 6.) compute DFTs of data_RF along first MATLAB dimension; choose appropriate order; in place
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// elapsed time for fft algorithm (ms)
	float time_dft = -1.0f;

	// print status
	if( verbosity > 1 ) mexPrintf(" %4s 2.) %s... ", "", "computing DFT");

	// place start event into the default stream
	checkCudaRTErrors( cudaEventRecord( start, 0 ) );

	// create plan for 1-D R2C DFT [batch size: N_signals]
	cufftHandle plan;
	checkCudaFFTErrors( cufftPlan1d( &plan, N_order_dft, CUFFT_R2C, N_signals ) );

	// execute DFT plan
	checkCudaFFTErrors( cufftExecR2C( plan, (cufftReal *) d_data_RF, (cufftComplex *) d_data_RF ) );

	// destroy plan after execution
	checkCudaFFTErrors( cufftDestroy( plan ) );

	// place stop event into the default stream
	checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
	checkCudaRTErrors( cudaEventSynchronize( stop ) );

	// compute elapsed time
	checkCudaRTErrors( cudaEventElapsedTime( &time_dft, start, stop ) );

	// print status
	if( verbosity > 1 )
	{
		mexPrintf( "done! (%.3f ms)\n", time_dft );
		// flush cache
		mexEvalString( "drawnow;" );
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 7.) copy addresses of __device__ window functions to host
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 11.) window function for the receive apodization
	t_window_ptr h_windows[ N_WINDOWS ];    // host-side function pointers to __device__ window functions
	checkCudaRTErrors( cudaMemcpyFromSymbol( &h_windows, d_windows, N_WINDOWS * sizeof( t_window_ptr ) ) );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 8.) compute image
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// elapsed time for kernel execution (ms)
	float time_kernel = -1.0f;

	// parallelization
	dim3 threadsPerBlock( N_THREADS_X, N_THREADS_Y );  // use max number of threads per block for additions
	dim3 numBlocks( N_blocks_x, N_blocks_z );

	// place start event into the default stream
	checkCudaRTErrors( cudaEventRecord( start, 0 ) );

	// check type of F-number
	if( f_number_tx_constant && f_number_rx_constant )
	{

		//--------------------------------------------------------------------------------------------------------------------------------------
		// a) fixed F-numbers
		//--------------------------------------------------------------------------------------------------------------------------------------
		// iterate tx elements
		for( int k_tx = 0; k_tx < N_el_tx; k_tx++ )
		{

			// compute propagation distances of incident spherical wave [ N_voxels x 1 ] -> move elsewhere and re-use?
			das_kernel_distances_sw<<<numBlocks, threadsPerBlock>>>( d_distances_tx,
															 		 d_positions_x, d_positions_z, N_positions_x, N_positions_z,
															 		 pos_rx_ctr_x[ k_tx ], 0.0f );

			// iterate rx elements
			for( int k_rx = 0; k_rx < N_el_rx; k_rx++ )
			{
				// compute DAS image for current rx channel
				// das_F_number_constant<<<numBlocks, threadsPerBlock>>>( d_image, d_weights_tx, d_weights_rx,
				// 											   d_distances_tx, argument_coefficient_dbl, argument_add_dbl,
				// 											   ( (cufftComplex *) d_data_RF ) + k_tx * N_el_rx * N_samples_dft + k_rx * N_samples_dft,
				// 											   d_positions_x, d_positions_z, N_positions_x, N_positions_z,
				// 											   pos_rx_ctr_x[ k_rx ], element_pitch_dbl, d_pos_rx_lb_x, d_pos_rx_ub_x, N_el_rx, M_el_rx,
				// 											   f_number_rx_values_dbl[ 0 ], h_windows[ window_rx_config.index ], window_rx_config.parameter,
				// 											   N_blocks_Omega_bp, index_f_lb, index_f_ub );
			} // for( int k_rx = 0; k_rx < N_el_rx; k_rx++ )

		} // for( int k_tx = 0; k_tx < N_el_tx; k_tx++ )

		// receive normalization
		if( normalization_rx )
		{
			// normalize image
			das_kernel_normalization<<<numBlocks, threadsPerBlock>>>( d_image, d_image, d_weights_rx, 0.0f, N_positions_x, N_positions_z );
		}

		// transmit normalization
		if( normalization_rx )
		{
			// normalize image
			das_kernel_normalization<<<numBlocks, threadsPerBlock>>>( d_image, d_image, d_weights_tx, 0.0f, N_positions_x, N_positions_z );
		}

	}
	else
	{

		//--------------------------------------------------------------------------------------------------------------------------------------
		// b) frequency-dependent F-numbers
		//--------------------------------------------------------------------------------------------------------------------------------------
		// allocate device memory for monofrequent DAS image
		cufftComplex* d_image_f = d_image;
		if( normalization_tx || normalization_rx ) checkCudaRTErrors( cudaMalloc( (void **) &d_image_f, size_bytes_DAS_image_RF_complex ) );

		// allocate device memory for matrix of weighted phase shifts
		cufftComplex* d_phase = NULL;
		checkCudaRTErrors( cudaMalloc( (void **) &d_phase, size_bytes_phase ) );

		// coefficients for cublasCgemm
		const t_float_complex_gpu gemm_alpha = make_cuFloatComplex( 1.0f, 0.0f );
		t_float_complex_gpu gemm_beta = make_cuFloatComplex( 1.0f, 0.0f );

		// create cuBLAS handle
		cublasHandle_t handle;
		checkCudaBLASErrors( cublasCreate( &handle ) );

		// iterate frequencies
		for( int index_f = index_f_lb; index_f <= index_f_ub; index_f++ )
		{

			// initialize monofrequent DAS image w/ zeros
			checkCudaRTErrors( cudaMemset( d_image_f, 0, size_bytes_DAS_image_RF_complex ) );

			// initialize weights w/ zeros
			checkCudaRTErrors( cudaMemset( d_weights_tx, 0, size_bytes_weights ) );

			// iterate tx elements
			for( int k_tx = 0; k_tx < N_el_tx; k_tx++ )
			{

				// compute propagation distances of incident wave [ N_voxels x 1 ] -> move elsewhere and re-use?
				das_kernel_distances_sw<<<numBlocks, threadsPerBlock>>>( d_distances_tx,
																 		 d_positions_x, d_positions_z, N_positions_x, N_positions_z,
																 		 pos_rx_ctr_x[ k_tx ], 0.0f );

				// compute matrix of complex-valued apodization weights [ N_voxels x N_elements ]
				saft_kernel_phase_shifts<<<numBlocks, threadsPerBlock>>>( d_phase, d_weights_tx, d_weights_rx,
																		  d_distances_tx, pos_rx_ctr_x[ k_tx ], index_f,
																		  d_positions_x, d_positions_z, N_positions_x, N_positions_z,
																		  d_pos_rx_ctr_x, element_pitch_dbl, d_pos_rx_lb_x, d_pos_rx_ub_x, N_el_rx, M_el_rx, N_blocks_rx,
																		  f_number_tx_values_dbl[ index_f - index_f_lb ], h_windows[ window_tx_config.index ], window_tx_config.parameter,
																		  f_number_rx_values_dbl[ index_f - index_f_lb ], h_windows[ window_rx_config.index ], window_rx_config.parameter,
																		  argument_coefficient_dbl, argument_add_dbl );

				// matrix-vector product
				// CUBLAS_OP_N: non-transpose operation / CUBLAS_OP_T: transpose operation / CUBLAS_OP_C: conjugate transpose operation
				// column-major format
				// CUBLAS_OP_T is required for d_data_RF because only a single frequency is read
				checkCudaBLASErrors(
					cublasCgemm( handle,
							CUBLAS_OP_N, CUBLAS_OP_T,
							N_positions_x * N_positions_z, 1, N_el_rx,
							&gemm_alpha, d_phase, N_positions_x * N_positions_z, ( (cufftComplex *) d_data_RF ) + k_tx * N_el_rx * N_samples_dft + index_f, N_samples_dft,
							&gemm_beta, d_image_f, N_positions_x * N_positions_z
						)
				);

			} // for( int k_tx = 0; k_tx < N_el_tx; k_tx++ )

			// receive normalization
			// TODO: use cuBLAS?
			if( normalization_rx )
			{
				// normalize the monofrequent image
				das_kernel_normalization<<<numBlocks, threadsPerBlock>>>( d_image_f, d_image_f, d_weights_rx, 0.0f, N_positions_x, N_positions_z );
			}

			// transmit normalization
			if( normalization_tx )
			{
				// normalize the monofrequent image and add to final image
				das_kernel_normalization<<<numBlocks, threadsPerBlock>>>( d_image, d_image_f, d_weights_tx, 1.0f, N_positions_x, N_positions_z );
			}

		} // for( int index_f = index_f_lb; index_f <= index_f_ub; index_f++ )

		// free device memory
		if( normalization_tx || normalization_rx ) checkCudaRTErrors( cudaFree( d_image_f ) );
		checkCudaRTErrors( cudaFree( d_phase ) );

		// destroy cuBLAS handle
		checkCudaBLASErrors( cublasDestroy( handle ) );

	} // if( f_number_tx_constant && f_number_rx_constant )

	// place stop event into the default stream
	checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
	checkCudaRTErrors( cudaEventSynchronize( stop ) );

	// compute elapsed time
	checkCudaRTErrors( cudaEventElapsedTime( &time_kernel, start, stop ) );

	// flush cache
	mexEvalString( "drawnow;" );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 9.) copy DAS image and weights to host
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 1.) allocate memory in the MATLAB workspace for the outputs
	//------------------------------------------------------------------------------------------------------------------------------------------
	// complex-valued DAS image
	plhs[ 0 ] = mxCreateNumericMatrix( N_positions_z, N_positions_x, mxSINGLE_CLASS, mxCOMPLEX );
	float* const image_DAS_real = (float*) mxGetData( plhs[ 0 ] );
	float* const image_DAS_imag = (float*) mxGetImagData( plhs[ 0 ] );

	// sum of apodization weights for each grid point (host)
	t_float_gpu* h_weights_tx = NULL; t_float_gpu* h_weights_rx = NULL;
	if( nlhs >= 2 )
	{
		plhs[ 1 ] = mxCreateNumericMatrix( N_positions_z, N_positions_x, mxSINGLE_CLASS, mxREAL );
		h_weights_tx = (t_float_gpu*) mxGetData( plhs[ 1 ] );
	}
	if( nlhs >= 3 )
	{
		plhs[ 2 ] = mxCreateNumericMatrix( N_positions_z, N_positions_x, mxSINGLE_CLASS, mxREAL );
		h_weights_rx = (t_float_gpu*) mxGetData( plhs[ 2 ] );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 2.) copy data from device to host
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) allocate memory for results
	cufftComplex* image_DAS_flt_cpx = (cufftComplex *) mxMalloc( size_bytes_DAS_image_RF_complex );

	// 	mexPrintf(" %4s 4.) %s... ", "", "transfering data to host");
	// place start event into the default stream
	checkCudaRTErrors( cudaEventRecord( start, 0 ) );

	// b) copy results to host
	checkCudaRTErrors( cudaMemcpy( image_DAS_flt_cpx, d_image, size_bytes_DAS_image_RF_complex, cudaMemcpyDeviceToHost ) );

	// c) copy transmit weights to host
	if( nlhs >= 2 )
	{
		checkCudaRTErrors( cudaMemcpy( h_weights_tx, d_weights_tx, size_bytes_weights, cudaMemcpyDeviceToHost ) );
	}

	// d) copy receive weights to host
	if( nlhs >= 3 )
	{
		checkCudaRTErrors( cudaMemcpy( h_weights_rx, d_weights_rx, size_bytes_weights, cudaMemcpyDeviceToHost ) );
	}

	// place stop event into the default stream
	checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
	checkCudaRTErrors( cudaEventSynchronize( stop ) );

	// compute elapsed time for memory transfer device to host (ms)
	float time_transfer_to_host = -1.0f;
	checkCudaRTErrors( cudaEventElapsedTime( &time_transfer_to_host, start, stop ) );

	// mexPrintf("done! (%.3f ms)\n", time_transfer_to_host);

	// d) extract complex-valued voxels
	for( int l_z = 0; l_z < N_positions_z; l_z++ )
	{
		for( int l_x = 0; l_x < N_positions_x; l_x++ )
		{
			// copy
			image_DAS_real[ l_x * N_positions_z + l_z ] = 2 * image_DAS_flt_cpx[ l_x * N_positions_z + l_z ].x / N_order_dft;
			image_DAS_imag[ l_x * N_positions_z + l_z ] = 2 * image_DAS_flt_cpx[ l_x * N_positions_z + l_z ].y / N_order_dft;

		} // for( int l_x = 0; l_x < N_positions_x; l_x++ )
	} // for( int l_z = 0; l_z < N_positions_z; l_z++ )

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 10.) print performance statistics
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// total elapsed time (ms)
	const float time_total = time_transfer_to_device + time_dft + time_kernel + time_transfer_to_host;

	if( verbosity > 0 )
	{
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf(" performance statistics:\n");
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf("%4s %-34s  %-24s  %-8s\n", "", "step", "time", "count");
		mexPrintf("%4s %-34s  %-24s  %-8s\n", "", "----------------------------------", "------------------------", "--------");
		mexPrintf("%4s %-34s: % 10.3f ms (%6.2f %%)  %8d\n", "", "memory transfer: host to device", time_transfer_to_device, time_transfer_to_device / time_total * 100, 1);
		mexPrintf("%4s %-34s: % 10.3f ms (%6.2f %%)  %8d\n", "", "DFT algorithm", time_dft, time_dft / time_total * 100, 1);
		mexPrintf("%4s %-34s: % 10.3f ms (%6.2f %%)  %8d\n", "", "kernel execution (per kernel)", time_kernel / N_el_rx, time_kernel / (N_el_rx * time_total) * 100, N_el_rx);
		mexPrintf("%4s %-34s: % 10.3f ms (%6.2f %%)  %8d\n", "", "kernel execution (total)", time_kernel, time_kernel / time_total * 100, 1);
		mexPrintf("%4s %-34s: % 10.3f ms (%6.2f %%)  %8d\n", "", "memory transfer: device to host", time_transfer_to_host, time_transfer_to_host / time_total * 100, 1);
		mexPrintf("%4s %-34s  %-24s  %-8s\n", "", "----------------------------------", "------------------------", "--------");
		mexPrintf("%4s %-34s: % 10.3f ms (%6.2f %%)\n", "", "total", time_total, 100.00f);
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 11.) clean up memory
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) device memory
	//------------------------------------------------------------------------------------------------------------------------------------------
	checkCudaRTErrors( cudaFree( d_positions_x ) );
	checkCudaRTErrors( cudaFree( d_positions_z ) );
	checkCudaRTErrors( cudaFree( d_data_RF ) );
	checkCudaRTErrors( cudaFree( d_pos_rx_ctr_x ) );
	checkCudaRTErrors( cudaFree( d_pos_rx_lb_x ) );
	checkCudaRTErrors( cudaFree( d_pos_rx_ub_x ) );
	checkCudaRTErrors( cudaFree( d_distances_tx ) );
	checkCudaRTErrors( cudaFree( d_image ) );
	checkCudaRTErrors( cudaFree( d_weights_tx ) );
	checkCudaRTErrors( cudaFree( d_weights_rx ) );

	checkCudaRTErrors( cudaEventDestroy( start ) );
	checkCudaRTErrors( cudaEventDestroy( stop ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// b) host memory
	//------------------------------------------------------------------------------------------------------------------------------------------
	// checkCudaErrors( cudaFreeHost( buffer ) );
	mxFree( positions_x );
	mxFree( positions_z );
	mxFree( data_RF );
	mxFree( pos_rx_ctr_x );
	mxFree( pos_rx_lb_x );
	mxFree( pos_rx_ub_x );
	mxFree( image_DAS_flt_cpx );

} // void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
