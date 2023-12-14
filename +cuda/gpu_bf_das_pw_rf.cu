/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % gpu_bf_das_pw_rf.cu
 %
 % Form ultrasound images by using
 % the delay-and-sum (DAS) algorithm in
 % the Fourier domain.
 %
 % The Fourier domain permits
 %  1.) independent receive focusing of different frequencies,
 %  2.) exact corrections of the arrival times independent of the sampling rate, and
 %  3.) usage of frequency-dependent apodization weights.
 %
 % A uniform linear transducer array is assumed.
 %
 % -------------------------------------------------------------------------
 % USAGE:
 % -------------------------------------------------------------------------
 % minimal:
 % image = gpu_bf_das_pw_rf( positions_x, positions_z, data_RF, f_s, steering_angle, element_width, element_pitch, c_0 );
 %
 % maximal:
 % [ image, F_number_values ] = gpu_bf_das_pw_rf( positions_x, positions_z, data_RF, f_s, steering_angle, element_width, element_pitch, c_0, f_bounds, index_t0, window, F_number, normalization, index_gpu );
 %
 % -------------------------------------------------------------------------
 % INPUTS:
 % -------------------------------------------------------------------------
 % REQUIRED
 %  00.) positions_x:           lateral voxel positions (m)
 %  01.) positions_z:           axial voxel positions (m)
 %  02.) data_RF:               RF data (2d array; 1st dimension: time, 2nd dimension: array element index)
 %  03.) f_s:                   sampling rate of the RF data (Hz)
 %  04.) steering_angle:        steering angle of the incident plane wave (rad) [ ← = pi/2, ↓ = 0, → = -pi/2 ]
 %  05.) element_width:         element width of the uniform linear array (m)
 %  06.) element_pitch:         element pitch of the uniform linear array (m)
 %  07.) c_0:                   speed of sound (m/s)
 %
 % OPTIONAL
 %  08.) f_bounds:              frequency bounds (Hz) (1d array; 1st element: lower frequency bound, 2nd element: upper frequency bound)
 %  09.) index_t0:              time index of the sample extracted from the focused RF signal (1)
 %  10.) window:                window function for receive apodization ( object of class windows.window )
 %  11.) F_number:              receive F-number (  object of class f_numbers.f_number )
 %  12.) normalization:         normalization of the complex-valued apodization weights ( object of class normalizations.normalization )
 %  13.) index_gpu:             index of the GPU to be used (1)
 %  14.) verbosity:             verbosity level of output
 %
 % -------------------------------------------------------------------------
 % OUTPUTS:
 % -------------------------------------------------------------------------
 %  00.) image:                 complex-valued DAS image
 %  01.) weights:               value of the F-number for each frequency
 %
 % -------------------------------------------------------------------------
 % ABOUT:
 % -------------------------------------------------------------------------
 % author: Martin Schiffner
 % date: 2010-12-17
 % modified: 2023-12-14
 % All rights reserved!
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

// TODO: compute order of DFT!
// __device__ t_float_gpu distance_pw( t_float_gpu pos_focus_x, t_float_gpu pos_focus_z, t_float_gpu pos_rx_ctr_x, t_float_gpu e_steering_x, t_float_gpu e_steering_z )
// {

// }

#include <cuda.h>
#include <cufft.h>
#include <math.h>
#include "gpu_bf_das_pw_rf.cuh"

// window functions
#include "gpu_bf_windows.cu"

// error handling
#include "gpu_bf_error_handling.cu"

#ifdef _CHAR16T
#define CHAR16_T
#endif
#include "mex.h"

#define N_THREADS_X 8
#define N_THREADS_Y 8
#define N_THREADS_PER_BLOCK 64

#define REVISION "1.2"
#define DATE "2023-12-14"
#define BYTES_PER_MEBIBYTE 1048576
#define BYTES_PER_KIBIBYTE 1024

// kernels
#include "das_pw_kernel_constant.cu"						// fixed F-number
#include "das_pw_kernel_frequency_dependent.cu"				// frequency-dependent F-number, boxcar normalization
#include "das_pw_kernel_frequency_dependent_window.cu"		// frequency-dependent F-number, window-based normalization

// workaround to check MATLAB classes
bool isa( const mxArray* object, const char* type )
{
	//Create LHS/RHS arrays for calling MATLAB
	mxArray* lhs[ 1 ];
	mxArray* rhs[ 2 ];
	//Return value
	bool retVal;
	//Populate Inputs to MATLAB isa
	rhs[ 0 ] = mxDuplicateArray( object );	// make deep copy of object to avoid problems w/ const modifier
	rhs[ 1 ] = mxCreateString( type );
	//Call the MATLAB isa function
	mexCallMATLAB( 1, lhs, 2, rhs, "isa" );
	//Extract result
	retVal = mxIsLogicalScalarTrue( lhs[ 0 ] );
	//Cleanup
	mxDestroyArray( lhs[ 0 ] );
	mxDestroyArray( rhs[ 0 ] );
	mxDestroyArray( rhs[ 1 ] );
	//Done
	return retVal;
}

// compute F-numbers
double* compute_values( const mxArray* f_number, mxArray* element_pitch_norm )
{
	mxArray* lhs[ 1 ];
	mxArray* rhs[ 2 ];
	rhs[ 0 ] = mxDuplicateArray( f_number );	// make deep copy of F-number to avoid problems w/ const modifier
	rhs[ 1 ] = element_pitch_norm;

	// call method compute_values in MATLAB
	mexCallMATLAB( 1, lhs, 2, rhs, "compute_values" );

	// clean up
	mxDestroyArray( rhs[ 0 ] );
	mxDestroyArray( rhs[ 1 ] );

	// get values
	return ( double * ) mxGetData( lhs[ 0 ] );
}

// map MATLAB window classes to ids
int get_window_id( const mxArray* window )
{

	if( mxIsClass( window, "windows.boxcar" ) )
	{
		return 0;
	}
	else if( mxIsClass( window, "windows.hann" ) )
	{
		return 1;
	}
	else if( mxIsClass( window, "windows.tukey" ) )
	{
		// ensure existence of property fraction_cosine
		if( mxGetProperty( window, 0, "fraction_cosine" ) != NULL )
		{
			// read property fraction_cosine
			if( *( ( double * ) mxGetData( mxGetProperty( window, 0, "fraction_cosine" ) ) ) == 0.2 )
			{
				return 2;
			}
		}
	}
	else if( mxIsClass( window, "windows.triangular" ) )
	{
		return 4;
	}

	return -1;

} // int get_window_id( const mxArray* window )

// convert object to string
int string( const mxArray* object, char* string, mwSize strlen )
{
	// if( !mxChar( object ) )
	// local variables
	mxArray* lhs_1;
	mxArray* lhs_2;
	mxArray* rhs_1;
	rhs_1 = mxDuplicateArray( object );	// make deep copy of F-number to avoid problems w/ const modifier

	// call method string in MATLAB
	mexCallMATLAB( 1, &lhs_1, 1, &rhs_1, "string" );

	// call method char in MATLAB
	mexCallMATLAB( 1, &lhs_2, 1, &lhs_1, "char" );

	// clean up
	mxDestroyArray( rhs_1 );
	mxDestroyArray( lhs_1 );

	// check for mxChar array
	// if( !mxIsChar( lhs_2 ) ) return 1;

	// mxChar array to C-style string
	return mxGetString( lhs_2, string, strlen );

}

  //Matlab's String class is encapsulated,
  //use Matlab call to convert it to char array
//   mxArray *string_class[1], *char_array[1];
//   string_class[0] = pr;
//   mexCallMATLAB(1, char_array, 1, string_class, "char");
  
  //Parse the char array to create an std::string
//   int buflen = mxGetN(char_array[0])*sizeof(mxChar)+1;
//   char* buf = new char[buflen];
//   mxGetString(char_array[0],buf,buflen);
//   data = std::string(buf);
//   delete buf;

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
	if( N_devices < 1 ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:NoCUDADevices", "Could not find any CUDA-enabled devices!" );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) check and assign inputs
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// check for correct number of input arguments
	if( nrhs < 8 ) mexErrMsgIdAndTxt("gpu_bf_das_pw_rf:nrhs", "At least 8 inputs are required.");
	if( nrhs > 15 ) mexErrMsgIdAndTxt("gpu_bf_das_pw_rf:nrhs", "At most 15 inputs are allowed.");

	// check for correct number of output arguments
	if( nlhs > 1 ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:nlhs", "Only one output is allowed.");

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 0.) lateral voxel positions (m)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with two dimensions
	if( !( mxIsDouble( prhs[ 0 ] ) && !mxIsComplex( prhs[ 0 ] ) && ( mxGetNumberOfDimensions( prhs[ 0 ] ) == 2 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidLateralPositions", "Lateral positions must be a real-valued matrix of type double!" );
	}

	// b) lateral positions
	const double* const positions_x_dbl = ( double * ) mxGetData( prhs[ 0 ] );

	// c) number of positions on x-axis
	const int N_pos_lat_x = mxGetN( prhs[ 0 ] );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 1.) axial voxel positions (m)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with two dimensions
	if( !( mxIsDouble( prhs[ 1 ] ) && !mxIsComplex( prhs[ 1 ] ) && ( mxGetNumberOfDimensions( prhs[ 1 ] ) == 2 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidAxialPositions", "Axial positions must be a real-valued matrix of type double!" );
	}

	// b) axial positions
	const double* const positions_z_dbl = ( double * ) mxGetData( prhs[ 1 ] );

	// c) number of positions on z-axis
	const int N_pos_lat_z = mxGetN( prhs[ 1 ] );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 2.) RF data (2d array; 1st dimension: time, 2nd dimension: array element index)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with two dimensions
	if( !( mxIsDouble( prhs[ 2 ] ) && !mxIsComplex( prhs[ 2 ] ) && ( mxGetNumberOfDimensions( prhs[ 2 ] ) == 2 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidRFData", "RF data must be a real-valued matrix of type double!" );
	}

	// b) create pointer to real data in RF data
	const double* const data_RF_dbl = ( double * ) mxGetData( prhs[ 2 ] );

	// c) read dimensions of two-dimensional RF data array
	const mwSize* const RF_data_dim = mxGetDimensions( prhs[ 2 ] );

	const int N_samples_t = RF_data_dim[ 0 ];  // number of samples per RF line in the time domain
	const int N_el_rx = RF_data_dim[ 1 ];      // number of array elements

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 3.) sampling rate of the RF data (Hz)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 3 ] ) && !mxIsComplex( prhs[ 3 ] ) && mxIsScalar( prhs[ 3 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSamplingRate", "Sampling rate must be a real-valued scalar of type double!" );
	}

	// b) sampling rate
	const double f_s_dbl = *( ( double * ) mxGetData( prhs[ 3 ] ) );

	// c) ensure positive sampling rate
	if( f_s_dbl <= 0 ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSamplingRate", "Sampling rate must be positive!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 4.) steering angle of the incident plane wave (rad) [ ← = pi/2, ↓ = 0, → = -pi/2 ]
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 4 ] ) && !mxIsComplex( prhs[ 4 ] ) && mxIsScalar( prhs[ 4 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSteeringAngle", "Steering angle must be a real-valued scalar of type double!" );
	}

	// b) steering angle of incident plane wave
	const double steering_angle_dbl = *( ( double * ) mxGetData( prhs[ 4 ] ) );

	// c) ensure valid steering angle between -90° and 90°
	if( ( steering_angle_dbl <= - M_PI / 2 ) || ( steering_angle_dbl >= M_PI / 2 ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSteeringAngle", "Steering angle must be larger than negative pi over two and smaller than pi over two!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 5.) element width of the uniform linear array (m)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 5 ] ) && !mxIsComplex( prhs[ 5 ] ) && mxIsScalar( prhs[ 5 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidElementWidth", "Element width must be a real-valued scalar of type double!" );
	}

	// b) element width
	const double element_width_dbl = *( ( double * ) mxGetData( prhs[ 5 ] ) );

	// c) ensure validity
	if( !( ( element_width_dbl > 0 ) && mxIsFinite( element_width_dbl ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidElementWidth", "Element width must be positive and finite!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 6.) element pitch of the uniform linear array (m)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 6 ] ) && !mxIsComplex( prhs[ 6 ] ) && mxIsScalar( prhs[ 6 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidElementPitch", "Element pitch must be a real-valued scalar of type double!" );
	}

	// b) element pitch
	const double element_pitch_dbl = *( ( double * ) mxGetData( prhs[ 6 ] ) );

	// c) ensure validity
	if( !( ( element_pitch_dbl > element_width_dbl ) && mxIsFinite( element_pitch_dbl ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidElementPitch", "Element pitch must be finite and larger than element width!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 7.) speed of sound (m/s)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) ensure real-valued mxDOUBLE_CLASS with one element
	if( !( mxIsDouble( prhs[ 7 ] ) && !mxIsComplex( prhs[ 7 ] ) && mxIsScalar( prhs[ 7 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSoundSpeed", "Speed of sound must be a real-valued scalar of type double!" );
	}

	// b) speed of sound
	const double c_0_dbl = *( ( double * ) mxGetData( prhs[ 7 ] ) );

	// c) ensure validity
	if( !( ( c_0_dbl > 0 ) && mxIsFinite( c_0_dbl ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSoundSpeed", "Average speed of sound must be positive and finite!" );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 8.) frequency bounds (Hz) (1d array; 1st element: lower frequency bound, 2nd element: upper frequency bound)
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default values
	double f_lb_dbl = DBL_MIN;		// lower bound: smallest double greater than zero
	double f_ub_dbl = f_s_dbl / 2; 	// upper bound: Nyquist frequency

	// b) check presence of nonempty argument
	if( nrhs >= 9 && !mxIsEmpty( prhs[ 8 ] ) )
	{
		// ensure real-valued mxDOUBLE_CLASS with two elements
		if( !( mxIsDouble( prhs[ 8 ] ) && !mxIsComplex( prhs[ 8 ] ) && ( mxGetNumberOfElements( prhs[ 8 ] ) == 2 ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyBounds", "Frequency bounds must be a real-valued vector with two entries of type double!" );
		}

		// frequency bounds
		const double* const f_bounds_dbl = ( double * ) mxGetData( prhs[ 8 ] );

		// extract frequency bounds
		f_lb_dbl = f_bounds_dbl[ 0 ];
		f_ub_dbl = f_bounds_dbl[ 1 ];
	}

	// c) ensure positive lower bound and valid upper bound
	if( f_lb_dbl <= 0 ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyBound", "Lower frequency bound must be positive!" );
	if( f_ub_dbl <= f_lb_dbl ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyBound", "Upper frequency bound must exceed the lower frequency bound!" );
	if( f_ub_dbl > f_s_dbl / 2 ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyBound", "Upper frequency bound must not exceed the Nyquist frequency!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 9.) time index of the sample extracted from the focused RF signal
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	double index_t0_dbl = 0.0;

	// b) check specified value
	if( nrhs >= 10 && !mxIsEmpty( prhs[ 9 ] ) )
	{
		// ensure real-valued mxDOUBLE_CLASS with one element
		if( !( mxIsDouble( prhs[ 9 ] ) && !mxIsComplex( prhs[ 9 ] ) && mxIsScalar( prhs[ 9 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidTimeShift", "Time index must be a real-valued scalar of type double!" );
		}

		// read value
		index_t0_dbl = *( ( double * ) mxGetData( prhs[ 9 ] ) );
	}

	// c) ensure nonnegative and finite time index
	if( !( ( index_t0_dbl >= 0 ) && mxIsFinite( index_t0_dbl ) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidTimeShift", "Time index must be nonnegative!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 10.) window function for the receive apodization ( object of class windows.window )
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	int index_window_rx = 0;

	// b) check specified value
	if( nrhs >= 11 && !mxIsEmpty( prhs[ 10 ] ) )
	{
		// ensure single object of class windows.window
		if( !( isa( prhs[ 10 ], "windows.window" ) && mxIsScalar( prhs[ 10 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidWindow", "Receive window must be a single object of the class windows.window!");
		}

		// map classes to window id
		index_window_rx = get_window_id( prhs[ 10 ] ); // index of apodization function
	}

	// c) check validity of window id
	if( !( ( index_window_rx >= 0 ) && ( index_window_rx < N_WINDOWS ) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidWindow", "Unknown window function!" );

	// d) name of window function
	char window_rx_name[ 128 ];
	if( string( prhs[ 10 ], window_rx_name, 128 * sizeof( char ) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidWindow", "Error reading window name!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 11.) F-number
	//------------------------------------------------------------------------------------------------------------------------------------------
//TODO: create default F-number

	// a) ensure single object of class f_numbers.f_number
	if( !( isa( prhs[ 11 ], "f_numbers.f_number" ) && mxIsScalar( prhs[ 11 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFNumber", "F-number must be a single object of the class f_numbers.f_number!");
	}

	// b) check if F-number is fixed
	int F_number_constant = 0;	// indicates whether F-number is fixed
	if( isa( prhs[ 11 ], "f_numbers.constant" ) )
	{
		F_number_constant = 1;
	}

	// c) name of F-number
	char f_number_name[ 128 ];
	if( string( prhs[ 11 ], f_number_name, 128 * sizeof(char) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFNumber", "Error reading F-number name!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 12.) normalization of the apodization weights
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	int index_window_f = 0;

	// b) ensure scalar object of class normalizations.normalization
	if( !( isa( prhs[ 12 ], "normalizations.normalization" ) && mxIsScalar( prhs[ 12 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidNormalization", "Normalization must be a single object of the class normalizations.normalization!");
	}

	// c) map classes to window id
	if( mxIsClass( prhs[ 12 ], "normalizations.windowed" ) )
	{
		index_window_f = get_window_id( mxGetProperty( prhs[ 12 ], 0, "window" ) );
	}

	// d) check validity of window id
	if( !( ( index_window_f >= 0 ) && ( index_window_f < N_WINDOWS ) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyWindow", "Unknown frequency window function!" );

	mexPrintf( "index_window_f = %d\n", index_window_f );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 13.) index of the GPU
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	int index_gpu = 0;

	// b) check specified value
	if( nrhs >= 14 && !mxIsEmpty( prhs[ 13 ] ) )
	{
		// ensure real-valued mxDOUBLE_CLASS with one element
		if( !( mxIsDouble( prhs[ 13 ] ) && !mxIsComplex( prhs[ 13 ] ) && mxIsScalar( prhs[ 13 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft:InvalidCUDADevice", "Device index must be a real-valued scalar of type double!" );
		}

		// read device index
		index_gpu = (int) *( ( double * ) mxGetData( prhs[ 13 ] ) );
	}

	// c) early exit if selected GPU does not exist
	if( ( index_gpu < 0 ) || ( index_gpu >= N_devices ) ) mexErrMsgIdAndTxt( "gpu_bf_saft:InvalidCUDADevice", "Selected device does not exist!" );
	// assertion:  index_gpu >= 0 && index_gpu < N_devices

	// d) select device
	checkCudaRTErrors( cudaSetDevice( index_gpu ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 14.) verbosity level of output
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) default value
	char verbosity = 1;

	// b) check specified value
	if( nrhs >= 15 && !mxIsEmpty( prhs[ 14 ] ) )
	{
		// ensure real-valued mxDOUBLE_CLASS with one element
		if( !( mxIsDouble( prhs[ 14 ] ) && !mxIsComplex( prhs[ 14 ] ) && mxIsScalar( prhs[ 14 ] ) ) )
		{
			mexErrMsgIdAndTxt( "gpu_bf_saft:InvalidVerbosityLevel", "Verbosity level must be a real-valued scalar of type double!" );
		}

		// read verbosity level
		verbosity = (char) *( ( double * ) mxGetData( prhs[ 14 ] ) );
	}

	// c) early exit if selected verbosity level does not exist
	if( ( verbosity < 0 ) || ( verbosity >= 4 ) ) mexErrMsgIdAndTxt( "gpu_bf_saft:InvalidVerbosityLevel", "Selected verbosity level does not exist!" );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 2.) compute dependent parameters and display status information
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) propagation direction of steered plane wave
	//------------------------------------------------------------------------------------------------------------------------------------------
	// x component of unit wave vector
	const double e_steering_x_dbl = -sin( steering_angle_dbl );

	// z component of unit wave vector
	const double e_steering_z_dbl = cos( steering_angle_dbl );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// b) array geometry
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
	// c) reference position
	//------------------------------------------------------------------------------------------------------------------------------------------
	// set reference position for the incident steered plane wave
	const double pos_tx_ctr_x_ref_dbl = ( e_steering_x_dbl >= 0 ) ? - M_el_rx * element_pitch_dbl : M_el_rx * element_pitch_dbl;

	//------------------------------------------------------------------------------------------------------------------------------------------
	// d) compute order of DFT; employ zero-padding to achieve odd order
	//------------------------------------------------------------------------------------------------------------------------------------------
// TODO: determine time duration of focused signal
	// tof_min, tof_max, t_lb, t_ub for all emissions

	// maximum relative time shift in the electronic focusing ( zero padding in DFT )
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
	// e) compute parameters of bandpass filter
	//------------------------------------------------------------------------------------------------------------------------------------------
	// lower frequency index in the RF signals
	const int index_f_lb = ceil( f_lb_dbl * N_order_dft / f_s_dbl );

	// upper frequency index in the RF signals
	const int index_f_ub = floor( f_ub_dbl * N_order_dft / f_s_dbl );

	// check validity of index_f_lb and index_f_ub
	if( ( index_f_lb > index_f_ub ) || ( index_f_lb <= 0 ) || ( index_f_ub >= N_samples_dft ) )
	{
		mexPrintf("error: invalid frequency boundaries f_lb and f_ub!");
		return;
	}
	// assertion: index_f_lb <= index_f_ub, index_f_lb > 0, index_f_ub < N_samples_dft

	// number of non-redundant samples of DFT of order N_order_dft after bandpass filter
	const int N_samples_dft_bp  = index_f_ub - index_f_lb + 1; // computed number of non-redundant samples of DFT of order N_order_dft after bandpass filter

	//------------------------------------------------------------------------------------------------------------------------------------------
	// f) compute frequency-dependent F-number [ in MATLAB ]
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
		const double* const f_number_values_dbl = compute_values( prhs[ 11 ], element_pitch_norm_matlab );
	// }

	//------------------------------------------------------------------------------------------------------------------------------------------
	// g) compute parallelization settings
	//------------------------------------------------------------------------------------------------------------------------------------------
	// number of blocks to process in parallel
	const int N_blocks_x = ceil( ( (double) N_pos_lat_x ) / N_THREADS_X );
	const int N_blocks_z = ceil( ( (double) N_pos_lat_z ) / N_THREADS_Y );

	// determine relevant range of frequency indices
	const int N_blocks_Omega_bp = ceil( ((double) N_samples_dft_bp) / N_THREADS_PER_BLOCK );  // number of frequency blocks

	//------------------------------------------------------------------------------------------------------------------------------------------
	// h) print global status information
	//------------------------------------------------------------------------------------------------------------------------------------------
	if( verbosity > 0 )
	{
		mexPrintf( " %s\n", "================================================================================" );
		mexPrintf( " GPU BF Toolbox: DAS, steered PW, v.%s (%s), Status Information\n", REVISION, DATE );
		mexPrintf( " (c) 2010-2023, M. F. Schiffner. All rights reserved.\n" );
		mexPrintf( " %s\n", "================================================================================" );
		mexPrintf( "  %-20s: % 5d x %-5d %5s", "lattice geometry", N_pos_lat_x, N_pos_lat_z, "" );
		mexPrintf( " %-22s: %-8d\n", "number of samples", N_samples_t );
		mexPrintf( "  %-20s: %-11d %7s", "number rx channels", N_el_rx, "" );
		mexPrintf( " %-22s: %5.1f MHz\n", "sampling rate", f_s_dbl / 1000000 );
		mexPrintf( "  %-20s: %7.2f m/s %7s", "average sound speed", c_0_dbl, "" );
		mexPrintf( " %-22s: %5.1f MHz (%d)\n", "lower frequency bound", f_lb_dbl / 1000000, index_f_lb );
		mexPrintf( "  %-20s: %6.2f° (%4.1f;%4.1f)", "angle of incidence", steering_angle_dbl * 180 / M_PI, e_steering_x_dbl, e_steering_z_dbl );
		mexPrintf( " %-22s: %5.1f MHz (%d)\n", "upper frequency bound", f_ub_dbl / 1000000, index_f_ub );
		mexPrintf( "  %-20s: %7.2f mm %8s", "reference position", pos_tx_ctr_x_ref_dbl * 1000, "" );
		mexPrintf( " %-22s: %-6d\n", "number of frequencies", N_samples_dft_bp );
		mexPrintf( "  %-20s: %-6.2f %12s", "additional shift", index_t0_dbl, "" );
		mexPrintf( " %-22s: %s (%d)\n", "window", window_rx_name, index_window_rx );
		mexPrintf( "  %-20s: %s\n", "F-number", f_number_name );
		mexPrintf( "  %-20s: %s\n", "normalization", "test" );

		//mexPrintf(" %-22s: %-6.2f um\n", "maximum delta_x", delta_x_max * 1e6);
		//mexPrintf(" %-22s: %-6.2f um\n", "maximum delta_z", delta_z_max * 1e6);
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 3.) infer and display properties of selected device
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
	// 4.) allocate and initialize device memory / convert data to float and copy to the device
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// size of conversion buffer: use maximum of
		// a) d_positions_x:
		//    sizeof( t_float_gpu ) * N_pos_lat_x
		// b) d_positions_z:
		//    sizeof( t_float_gpu ) * N_pos_lat_z;
		// c) d_data_RF:
		//    sizeof( cufftComplex ) * N_samples_dft * N_el_rx
		// d) d_pos_rx_ctr_x, d_pos_rx_lb_x, d_pos_rx_ub_x:
		//    sizeof( t_float_gpu ) * N_el_rx
		// e) d_f_number_values:
		//    sizeof( float ) * N_samples_dft_bp
		// f) d_image:
		//    sizeof( cufftComplex ) * N_pos_lat_x * N_pos_lat_z
		// g) d_weights:
		//    sizeof( t_float_gpu ) * N_pos_lat_x * N_pos_lat_z
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

	// 	// place start event into the default stream
	// 	checkCudaRTErrors( cudaEventRecord( start, 0 ) );
	
	// 	// place stop event into the default stream
	// 	checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
	// 	checkCudaRTErrors( cudaEventSynchronize( stop ) );

	// 	// compute elapsed time
	float time_transfer_to_device = -1.0f;  // elapsed time for memory transfer host to device (ms)
	// 	checkCudaRTErrors( cudaEventElapsedTime( &time_transfer_to_device, start, stop ) );

	// 	mexPrintf("done! (%.3f ms)\n", time_transfer_to_device);

	// 	// flush cache
	// 	mexEvalString( "drawnow;" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 0.) lateral voxel positions
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes (FLOAT)
	const size_t size_bytes_positions_x = sizeof( t_float_gpu ) * N_pos_lat_x;

	// b) allocate memory
	t_float_gpu* positions_x = (t_float_gpu *) mxMalloc( size_bytes_positions_x );

	// c) convert lateral image coordinates
	for( int l_x = 0; l_x < N_pos_lat_x; l_x++ )
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
	const size_t size_bytes_positions_z = sizeof( t_float_gpu ) * N_pos_lat_z;

	// b) allocate memory
	t_float_gpu* positions_z = (t_float_gpu *) mxMalloc( size_bytes_positions_z );

	// c) conversion to selected data type
	for( int l_z = 0; l_z < N_pos_lat_z; l_z++ )
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
	const size_t size_bytes_RF_data = sizeof( cufftComplex ) * N_samples_dft * N_el_rx;
	const size_t size_bytes_RF_data_host = sizeof( t_float_gpu ) * N_samples_t * N_el_rx;

	// b) allocate memory
	t_float_gpu* data_RF = (t_float_gpu *) mxMalloc( size_bytes_RF_data_host );

	// c) conversion to selected data type
	for( int k_rx = 0; k_rx < N_el_rx; k_rx++ )
	{
		for( int l_z = 0; l_z < N_samples_t; l_z++ )
		{
			data_RF[ k_rx * N_samples_t + l_z ] = (t_float_gpu) data_RF_dbl[ k_rx * N_samples_t + l_z ];      
		}
	}

	// d) allocate device memory for RF data and its DFT [additional memory for in-place DFT]
	t_float_gpu* d_data_RF = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_data_RF, size_bytes_RF_data ) );

	// e) initialize RF data (necessary because of zero-padding)
	checkCudaRTErrors( cudaMemset( d_data_RF, 0, size_bytes_RF_data ) );

	// f) copy data [correctly pad RF data to allow in place transform]
	const size_t size_bytes_signal  = sizeof( t_float_gpu ) * N_samples_t;
	for( int k_rx = 0; k_rx < N_el_rx; k_rx++ )
	{
		checkCudaRTErrors( cudaMemcpy( d_data_RF + k_rx * N_samples_dft * 2, &( data_RF[ k_rx * N_samples_t ] ), size_bytes_signal, cudaMemcpyHostToDevice ) );
	}
	// N_order_dft / 2 + 1 | data_RF + k_rx * N_samples_t

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 3.) array geometry
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes
	const size_t size_bytes_positions_rx = sizeof( t_float_gpu ) * N_el_rx;

	// b) allocate device memory for RF data and its DFT [additional memory for in-place DFT]
	t_float_gpu* d_pos_rx_ctr_x = NULL; t_float_gpu* d_pos_rx_lb_x = NULL; t_float_gpu* d_pos_rx_ub_x = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_ctr_x, size_bytes_positions_rx) );
	checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_lb_x, size_bytes_positions_rx) );
	checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_ub_x, size_bytes_positions_rx) );

	// c) copy data
	checkCudaRTErrors( cudaMemcpy( d_pos_rx_ctr_x, pos_rx_ctr_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );
	checkCudaRTErrors( cudaMemcpy( d_pos_rx_lb_x, pos_rx_lb_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );
	checkCudaRTErrors( cudaMemcpy( d_pos_rx_ub_x, pos_rx_ub_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 4.) F-number
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes (FLOAT)
	const size_t size_bytes_f_numbers = sizeof( t_float_gpu ) * N_samples_dft_bp;

	// b) allocate memory
	t_float_gpu* f_number_values = (t_float_gpu*) mxMalloc( size_bytes_f_numbers );

	// c) convert frequency-dependent F-number
	for( int index_f = 0; index_f < N_samples_dft_bp; index_f++ )
	{
		f_number_values[ index_f ] = (t_float_gpu) f_number_values_dbl[ index_f ];
		// mexPrintf( "f_number_values[ %d ] = %.2f\n", index_f, f_number_values[ index_f ] );
	}

	// d) allocate device memory
	t_float_gpu* d_f_number_values = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_f_number_values, size_bytes_f_numbers ) );

	// e) copy data
	checkCudaRTErrors( cudaMemcpy( d_f_number_values, f_number_values, size_bytes_f_numbers, cudaMemcpyHostToDevice ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 5.) TODO: arguments
	//------------------------------------------------------------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 6.) image
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes
	const size_t size_bytes_DAS_image_RF_complex = sizeof( cufftComplex ) * N_pos_lat_x * N_pos_lat_z;

	// b) allocate device memory
	cufftComplex* d_image = NULL;
	checkCudaRTErrors( cudaMalloc( (void **) &d_image, size_bytes_DAS_image_RF_complex ) );

	// c) initialize complex-valued DAS image
	checkCudaRTErrors( cudaMemset( d_image, 0, size_bytes_DAS_image_RF_complex ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 7.) weights
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) compute required size for array in bytes
	size_t size_bytes_weights = sizeof( t_float_gpu ) * N_pos_lat_x * N_pos_lat_z;

	// b) allocate device memory
	t_float_gpu* d_weights = NULL;
	checkCudaRTErrors( cudaMalloc( (void**) &d_weights, size_bytes_weights ) );

	// c) initialize weights w/ zeros
	checkCudaRTErrors( cudaMemset( d_weights, 0, size_bytes_weights ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 7.) memory status
	//------------------------------------------------------------------------------------------------------------------------------------------
	// total memory consumption
	const size_t size_bytes_total = size_bytes_positions_x + size_bytes_positions_z + size_bytes_RF_data + 3 * size_bytes_positions_rx + size_bytes_f_numbers + size_bytes_DAS_image_RF_complex + size_bytes_weights;

	// get available amount of global memory on GPU device
	size_t mem_available = 0;
	checkCudaRTErrors( cudaMemGetInfo( &mem_available, &(deviceProp.totalGlobalMem) ) );

	if( verbosity > 0 )
	{
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf(" memory consumption (available global memory: %-6.2f MiByte (%6.2f %%))\n", ((double) mem_available) / BYTES_PER_MEBIBYTE, ((double) mem_available) * 100 / deviceProp.totalGlobalMem);
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf(" %4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "RF signals", ((double) size_bytes_RF_data) / BYTES_PER_KIBIBYTE, ((double) size_bytes_RF_data) * 100 / size_bytes_total);
		mexPrintf(" %4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "lateral grid positions", ((double) size_bytes_positions_x) / BYTES_PER_KIBIBYTE, ((double) size_bytes_positions_x) * 100 / size_bytes_total);
		mexPrintf(" %4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "axial grid positions", ((double) size_bytes_positions_z) / BYTES_PER_KIBIBYTE, ((double) size_bytes_positions_z) * 100 / size_bytes_total);
		mexPrintf(" %4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "array geometry", ((double) 3.0 * size_bytes_positions_rx) / BYTES_PER_KIBIBYTE, ((double) 3.0 * size_bytes_positions_rx) * 100 / size_bytes_total);
		mexPrintf(" %4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "F-numbers", ((double) size_bytes_f_numbers) / BYTES_PER_KIBIBYTE, ((double) size_bytes_f_numbers) * 100 / size_bytes_total);
		mexPrintf(" %4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "DAS image", ((double) size_bytes_DAS_image_RF_complex) / BYTES_PER_KIBIBYTE, ((double) size_bytes_DAS_image_RF_complex) * 100 / size_bytes_total);        
		mexPrintf(" %4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "weights", ((double) size_bytes_weights) / BYTES_PER_KIBIBYTE, ((double) size_bytes_weights) * 100 / size_bytes_total);
		mexPrintf(" %4s %s  %s\n", "", "-----------------------------", "----------------------------");
		mexPrintf(" %4s %-29s: % 10.2f KiByte (%6.2f %%)\n", "", "total", ((double) size_bytes_total) / BYTES_PER_KIBIBYTE, 100.0);
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 5.) compute DFTs of data_RF along first MATLAB dimension; choose appropriate order; in place
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// elapsed time for fft algorithm (ms)
	float time_dft = -1.0f;

	// print status
	if( verbosity > 1 ) mexPrintf(" %4s 2.) %s... ", "", "computing DFT");

	// create events for time-measurement
	cudaEvent_t start, stop;
	checkCudaRTErrors( cudaEventCreate( &start ) );
	checkCudaRTErrors( cudaEventCreate( &stop ) );

	// place start event into the default stream
	checkCudaRTErrors( cudaEventRecord( start, 0 ) );

	// create plan for 1-D R2C DFT [batch size: N_el_rx]
	cufftHandle plan;
	checkCudaFFTErrors( cufftPlan1d( &plan, N_order_dft, CUFFT_R2C, N_el_rx ) );

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
	// 6.) copy addresses of __device__ window functions to host
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 11.) window function for the receive apodization
	t_window_ptr h_windows[ N_WINDOWS ];    // host-side function pointers to __device__ window functions
	checkCudaRTErrors( cudaMemcpyFromSymbol( &h_windows, d_windows, N_WINDOWS * sizeof( t_window_ptr ) ) );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 7.) compute DAS image on GPU
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 	mexPrintf(" %4s 3.) %s... ", "", "computing DAS image on GPU");
	float time_kernel = -1.0f;              // elapsed time for kernel execution (ms)
	dim3 threadsPerBlock(N_THREADS_X, N_THREADS_Y);  /* use max number of threads per block for additions*/
	dim3 numBlocks(N_blocks_x, N_blocks_z);

	// place start event into the default stream
	checkCudaRTErrors( cudaEventRecord( start, 0 ) );

	// TODO: kernel to compute TOFs; d_arguments: N_el_rx * N_positions_x * N_positions_z

	// iterate rx channels
	for( int k_rx = 0; k_rx < N_el_rx; k_rx++ )
	{

		// compute DAS image for current tx/rx configuration
		if( F_number_constant )
		{
			das_F_number_constant<<<numBlocks, threadsPerBlock>>>( d_image, d_weights,
																   ((cufftComplex *) d_data_RF) + k_rx * N_samples_dft,
																   d_positions_x, d_positions_z, pos_rx_ctr_x[ k_rx ], d_pos_rx_lb_x, d_pos_rx_ub_x,
																   f_number_values[ 0 ],
																   h_windows[ index_window_rx ],
																   argument_coefficient_dbl, N_blocks_Omega_bp, index_f_lb, index_f_ub,
																   N_pos_lat_x, N_pos_lat_z, e_steering_x_dbl, e_steering_z_dbl, pos_tx_ctr_x_ref_dbl,
																   argument_add_dbl,
																   element_pitch_dbl, M_el_rx, N_el_rx );
		}
		else
		{
			if( index_window_f == 0)
			{
				// boxcar window along bandwidth
				das_F_number_frequency_dependent<<<numBlocks, threadsPerBlock>>>( d_image,
																				  ((cufftComplex *) d_data_RF) + k_rx * N_samples_dft,
																				  d_positions_x, d_positions_z,
																				  d_pos_rx_ctr_x, d_pos_rx_lb_x, d_pos_rx_ub_x,
																				  d_f_number_values,
																				  h_windows[ index_window_rx ],
																				  argument_coefficient_dbl, N_blocks_Omega_bp, index_f_lb, index_f_ub,
																				  N_pos_lat_x, N_pos_lat_z, e_steering_x_dbl, e_steering_z_dbl, pos_tx_ctr_x_ref_dbl,
																				  argument_add_dbl, k_rx,
																				  element_pitch_dbl, M_el_rx, N_el_rx );
			}
			else
			{
				// // other window
				// das_F_number_frequency_dependent_window<<<numBlocks, threadsPerBlock>>>( d_image,
				// 																		 ((cufftComplex *) d_data_RF) + k_rx * N_samples_dft,
				// 																		 d_positions_x, d_positions_z,
				// 																		 d_pos_rx_ctr_x, d_pos_rx_lb_x, d_pos_rx_ub_x,
				// 																		 d_f_number_values,
				// 																		 h_windows[ index_window_rx ], h_windows[ index_window_f ],
				// 																		 argument_coefficient_dbl, N_blocks_Omega_bp, index_f_lb, N_samples_dft_bp,
				// 																		 N_pos_lat_x, N_pos_lat_z, e_steering_x, e_steering_z, pos_tx_ctr_x_ref_dbl,
				// 																		 argument_add_dbl, f_s_over_c_0_dbl, k_rx,
				// 																		 ( t_float_gpu ) element_pitch_dbl, ( t_float_gpu ) M_el_rx, N_el_rx );
			} // if( index_window_f == 0)
		} // if( F_number_constant )

	} // for( k_rx = 0; k_rx < N_el_rx; k_rx++ )

	// place stop event into the default stream
	checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
	checkCudaRTErrors( cudaEventSynchronize( stop ) );

	// compute elapsed time
	checkCudaRTErrors( cudaEventElapsedTime( &time_kernel, start, stop ) );

	// mexPrintf("done! (%.3f ms)\n", time_kernel);

	// flush cache
	mexEvalString( "drawnow;" );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 8.) copy DAS image and weights to host
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 1.) allocate memory in the MATLAB workspace for the output
	//------------------------------------------------------------------------------------------------------------------------------------------
	// complex-valued DAS image
	plhs[ 0 ] = mxCreateNumericMatrix( N_pos_lat_z, N_pos_lat_x, mxSINGLE_CLASS, mxCOMPLEX );
	float* const image_DAS_real = (float*) mxGetData( plhs[ 0 ] );
	float* const image_DAS_imag = (float*) mxGetImagData( plhs[ 0 ] );

	// sum of apodization weights for each grid point (host)
	plhs[ 1 ] = mxCreateNumericMatrix( N_pos_lat_z, N_pos_lat_x, mxSINGLE_CLASS, mxREAL );
	t_float_gpu* const h_weights = (t_float_gpu*) mxGetData( plhs[ 1 ] );

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

	// c) copy weights to host
	if( F_number_constant )
	{
		checkCudaRTErrors( cudaMemcpy( h_weights, d_weights, size_bytes_weights, cudaMemcpyDeviceToHost ) );
	}

	// place stop event into the default stream
	checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
	checkCudaRTErrors( cudaEventSynchronize( stop ) );

	// compute elapsed time for memory transfer device to host (ms)
	float time_transfer_to_host = -1.0f;
	checkCudaRTErrors( cudaEventElapsedTime( &time_transfer_to_host, start, stop ) );

	// 	mexPrintf("done! (%.3f ms)\n", time_transfer_to_host);

	// d) extract complex-valued pixels and divide by weights (if appropriate)
	for( int l_z = 0; l_z < N_pos_lat_z; l_z++ )
	{
		for( int l_x = 0; l_x < N_pos_lat_x; l_x++ )
		{

			//
			image_DAS_real[ l_x * N_pos_lat_z + l_z ] = image_DAS_flt_cpx[ l_x * N_pos_lat_z + l_z ].x;
			image_DAS_imag[ l_x * N_pos_lat_z + l_z ] = image_DAS_flt_cpx[ l_x * N_pos_lat_z + l_z ].y;

			// normalize pixels for constant F-number
			if( F_number_constant )
			{
				// check how many rx channels contributed to current pixel
				if( h_weights[ l_x * N_pos_lat_z + l_z ] > 0 )
				{
					// image_DAS_real[ l_x * N_pos_lat_z + l_z ] = image_DAS_real[ l_x * N_pos_lat_z + l_z ] / ( N_order_dft * h_weights[ l_x * N_pos_lat_z + l_z ] );
					// image_DAS_imag[ l_x * N_pos_lat_z + l_z ] = image_DAS_imag[ l_x * N_pos_lat_z + l_z ] / ( N_order_dft * h_weights[ l_x * N_pos_lat_z + l_z ] );
					image_DAS_real[ l_x * N_pos_lat_z + l_z ] = image_DAS_real[ l_x * N_pos_lat_z + l_z ] / ( h_weights[ l_x * N_pos_lat_z + l_z ] );
					image_DAS_imag[ l_x * N_pos_lat_z + l_z ] = image_DAS_imag[ l_x * N_pos_lat_z + l_z ] / ( h_weights[ l_x * N_pos_lat_z + l_z ] );
				}
				// else
				// {
				// 	// no rx channel contributed to current pixel
				// 	image_DAS_real[ l_x * N_pos_lat_z + l_z ] = 0;
				// 	image_DAS_imag[ l_x * N_pos_lat_z + l_z ] = 0;
				// }
			} // if( F_number_constant )

		} // for( int l_x = 0; l_x < N_pos_lat_x; l_x++ )
	} // for( int l_z = 0; l_z < N_pos_lat_z; l_z++ )

	// mxDestroyArray( lhs[ 0 ] );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 9.) print performance statistics
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
	// 10.) clean up memory
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
	checkCudaRTErrors( cudaFree( d_f_number_values ) );
	checkCudaRTErrors( cudaFree( d_image ) );
	checkCudaRTErrors( cudaFree( d_weights ) );

	checkCudaRTErrors( cudaEventDestroy( start ) );
	checkCudaRTErrors( cudaEventDestroy( stop ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// b) host memory
	//------------------------------------------------------------------------------------------------------------------------------------------
	// checkCudaErrors( cudaFreeHost( buffer ) );
	mxFree( f_number_values );
	mxFree( data_RF );
	mxFree( positions_z );
	mxFree( positions_x );
	mxFree( pos_rx_ub_x  );
	mxFree( pos_rx_lb_x  );
	mxFree( pos_rx_ctr_x );
	mxFree( image_DAS_flt_cpx );

} // void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
