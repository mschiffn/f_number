/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % gpu_bf_das_pw_rf.cu
 %
 % delay-and-sum (DAS) beamforming in the temporal Fourier domain using
 % the RF signals obtained from the excitation by
 % a steered plane wave
 % 
 % The calling syntax is:
 %
 % image_DAS
 % = gpu_bf_das_pw_rf
 % (
 % 00.) pos_lat_x,				lateral positions of the image voxels (m)
 % 01.) pos_lat_z,				axial positions of the image voxels (m)
 % 02.) data_RF,				real-valued RF signals (2d array; 1st dimension: time, 2nd dimension: array element index)
 % 03.) f_s,					sampling rate of RF signals (Hz)
 % 04.) theta_incident,			steering angle of the incident plane wave (rad)
 % 05.) element_width,			element width of the linear array (m)
 % 06.) element_pitch,			element pitch of the linear array (m)
 % 07.) f_lb,					lower frequency bound in the RF signals (Hz)
 % 08.) f_ub,					upper frequency bound in the RF signals (Hz)
 % 09.) c_0,					average small-signal sound speed (m/s)
 % 10.) N_samples_shift_add,	number of additional samples to shift (1)
 % 11.) window,					window function for the receive apodization weights ( object of class windows.window )
 % 12.) F_number,				receive F-number ( object of class f_numbers.f_number )
 % 13.) normalization,			normalization of the receive apodization weights ( object of class normalizations.normalization )
 % 14.) index_gpu				index of the GPU to be used (1)
 % )
 %
 % author: Martin Schiffner
 % date: 2010-12-17
 % modified: 2022-03-03
 % All rights reserved!
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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

#define REVISION "1.1"
#define DATE "2022-02-18"
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
	// 1.) define variables on host
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 0.) pointer to one-dimensional arrays with rx focus locations on x-axis
	double* pos_lat_x_dbl = NULL;		// double precision
	float* pos_lat_x = NULL;			// single precision

	// 1.) pointer to one-dimensional arrays with rx focus locations on z-axis
	double* pos_lat_z_dbl = NULL;		// double precision
	float* pos_lat_z = NULL;			// single precision

	// 2.) pointer to two-dimensional array of real-valued RF signals
	double *data_RF_dbl = NULL;			// double precision
	float *data_RF = NULL;				// single precision

	// 3.) temporal sampling frequency
	double f_s_dbl = 0.0;				// double precision

	// 4.) angle of incident plane wave
	double theta_incident_dbl = 0.0;	// double precision

	// x component of unit wave vector
	double cos_theta_incident_dbl = 0.0;	// double precision
	float cos_theta_incident_flt = 0.0f;	// single precision

	// z component of unit wave vector
	double sin_theta_incident_dbl = 0.0;	// double precision
	float sin_theta_incident_flt = 0.0f;	// single precision

	// 5.) width of array elements
	double element_width_dbl = 0.0;		// double precision

	// 6.) width of array elements
	double element_pitch_dbl = 0.0;		// double precision

	// 7.) lower bound on the temporal frequency
	double f_lb_dbl = 0.0;				// double precision
	int index_f_lb = 0;

	// 8.) upper bound on the temporal frequency
	double f_ub_dbl = 0.0;				// double precision
	int index_f_ub = 0;

	// 9.) average small-signal sound speed
	double c_0_dbl = 0.0;				// double precision

	// 10.) number of samples added to computed time delay
	double N_samples_shift_add_dbl = 0.0;	// double precision (default value: 0)
	float N_samples_shift_add_flt = 0.0f;	// single precision

	// 11.) window function for the receive apodization
	t_window_ptr h_windows[ N_WINDOWS ];	// host-side function pointers to __device__ window functions
	int index_window_rx = 0;				// index of apodization function (default value: 0)
	char window_rx_name[ 128 ];				// name of window function

	// 12.) frequency-dependent F-number
	double* f_number_values_dbl = NULL;		// double precision
	float* f_number_values_flt = NULL;		// single precision
	char f_number_name[ 128 ];				// name of F-number

	// 13.) index of the frequency window function
	int index_window_f = 0;				// default value

	// 14.) index of GPU to be used
	int index_gpu = 0;					// default value

	// pointers to one-dimensional array with rx locations on x-axis
	float* pos_rx_ctr_x = NULL;			// single precision
	float* pos_rx_lb_x = NULL;			// single precision
	float* pos_rx_ub_x = NULL;			// single precision

	// x position of reference element (tx at t = 0)
	double pos_tx_x_ref_dbl = 0.0;		// double precision
	float pos_tx_ctr_x_ref_flt = 0.0f;		// single precision

	// frequency axis
	double* element_pitch_norm_dbl = NULL;

	// resulting DAS image
	cufftComplex* image_DAS_flt_cpx;	// complex data, single precision
	float* image_DAS_real;              // real data, single precision
	float* image_DAS_imag;              // imaginary data, single precision

	// sum of apodization weights for each grid point
	float* h_weights = NULL;			// sum of apodization weights for each grid point (host)

	// quotients
	float f_s_over_c_0_flt;				// single precision

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) define variables on device
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	cufftComplex* d_image_flt_cpx;		// pointer to DAS image on GPU, complex data, single precision
	float* data_RF_flt_gpu;				// pointer to RF data on GPU, single precision
	float* d_pos_lat_x_flt;				// pointer to one-dimensional arrays with rx focus locations on x-axis
	float* d_pos_lat_z_flt;				// pointer to one-dimensional arrays with rx focus locations on z-axis

	t_float_gpu* d_pos_rx_ctr_x = NULL;		// array elements
	t_float_gpu* d_pos_rx_lb_x = NULL;
	t_float_gpu* d_pos_rx_ub_x = NULL;

	t_float_gpu* d_weights = NULL;		// pointer to weights for image pixels on GPU

	t_float_gpu* d_f_number_values = NULL;		//

    // timing variables for performance statistics
    cudaEvent_t start, stop;
    float time_transfer_to_device = -1.0f;  // elapsed time for memory transfer host to device (ms)
    float time_dft = -1.0f;                 // elapsed time for fft algorithm (ms)
    float time_kernel = -1.0f;              // elapsed time for kernel execution (ms)
    float time_transfer_to_host = -1.0f;    // elapsed time for memory transfer device to host (ms)
    float time_total = 0.0f;                // elapsed time total (ms)
    
    double argument_factor_dbl;
    float argument_factor_flt;

	const mwSize *RF_data_dim;  // dimensions of two-dimensional RF data array
	int N_el_rx = 0;			// number of receive elements
	double M_el_rx = 0;			// half-number of receive elements

	int N_pos_lat_x = 0;        // number of positions on x-axis
	int N_pos_lat_z = 0;        // number of positions on z-axis

	int N_samples_t = 0;        // number of samples per RF line

	int N_order_dft = 0;        // order of DFT

	int N_samples_dft = 0;      // computed number of non-redundant samples of DFT of order N_order_dft

	int N_samples_dft_bp = 0;   // computed number of non-redundant samples of DFT of order N_order_dft after bandpass filter

	int N_blocks_Omega_bp = 0;	// number of frequency blocks

	size_t size_bytes_DAS_image_RF_complex, size_bytes_RF_data, size_bytes_positions_x, size_bytes_positions_z, size_bytes_positions_rx, size_bytes_f_numbers, size_bytes_weights, size_bytes_signal, mem_available, size_bytes_total;

	int cufft_batch_size = 10;  // batch size for 1-D R2C transform in cuFFT library

	int l_x = 0;                // loop variables
	int l_z = 0;               
	int k_rx = 0;

	int N_blocks_x = 0;
	int N_blocks_z = 0;
	int F_number_constant = 0;	// indicates whether F-number is fixed

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 2.) check and assign input arguments
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// check for correct number of input arguments
	if( nrhs < 10 ) mexErrMsgIdAndTxt("gpu_bf_das_pw_rf:nrhs", "At least 10 inputs are required.");
	if( nrhs > 15 ) mexErrMsgIdAndTxt("gpu_bf_das_pw_rf:nrhs", "At most 15 inputs are allowed.");

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 0.) lateral positions
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 0 ] ) && ( mxGetNumberOfDimensions( prhs[ 0 ] ) == 2 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidLateralPositions", "Lateral positions must be two-dimensional and double!" );
	}

	// b) lateral positions
	pos_lat_x_dbl = ( double * ) mxGetData( prhs[ 0 ] );
	N_pos_lat_x = mxGetN( prhs[ 0 ] );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 1.) axial positions
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 1 ] ) && ( mxGetNumberOfDimensions( prhs[ 1 ] ) == 2 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidAxialPositions", "Axial positions must be two-dimensional and double!" );
	}

	// b) axial positions
	pos_lat_z_dbl = ( double * ) mxGetData( prhs[ 1 ] );
	N_pos_lat_z = mxGetN( prhs[ 1 ] );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 2.) RF data
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 2 ] ) && ( mxGetNumberOfDimensions( prhs[ 2 ] ) == 2 ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidRFData", "RF data must be double matrix!" );
	}

	// b) create pointer to real data in RF data
	data_RF_dbl = ( double * ) mxGetData( prhs[ 2 ] );

	// c) read dimensions of RF data
	RF_data_dim = mxGetDimensions( prhs[ 2 ] );

	N_samples_t = RF_data_dim[ 0 ];		// number of samples in the time domain
	N_el_rx     = RF_data_dim[ 1 ];		// number of array elements
	M_el_rx = ( N_el_rx - 1.0 ) / 2.0;	// half-number of receive elements

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 3.) sampling rate
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 3 ] ) && mxIsScalar( prhs[ 3 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSamplingRate", "Sampling rate must be double and scalar!" );
	}

	// b) sampling rate
	f_s_dbl = *( ( double * ) mxGetData( prhs[ 3 ] ) );

	// c) ensure positive sampling rate
	if( f_s_dbl <= 0 ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSamplingRate", "Sampling rate must be positive!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 4.) steering angle of incident plane wave
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 4 ] ) && mxIsScalar( prhs[ 4 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSteeringAngle", "Steering angle must be double and scalar!" );
	}

	// b) steering angle of incident plane wave
	theta_incident_dbl = *( ( double * ) mxGetData( prhs[ 4 ] ) );

	// c) ensure valid steering angle
	if( ( theta_incident_dbl <= 0 ) || ( theta_incident_dbl >= M_PI ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSteeringAngle", "Steering angle must be positive and smaller than pi!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 5.) element width
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 5 ] ) && mxIsScalar( prhs[ 5 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidElementWidth", "Element width must be double and scalar!" );
	}

	// b) element width
	element_width_dbl = *( ( double * ) mxGetData( prhs[ 5 ] ) );

	// c) ensure valid element width
	if( element_width_dbl <= 0 ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidElementWidth", "Element width must be positive!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 6.) element pitch
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 6 ] ) && mxIsScalar( prhs[ 6 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidElementPitch", "Element pitch must be double and scalar!" );
	}

	// b) element pitch
	element_pitch_dbl = *( ( double * ) mxGetData( prhs[ 6 ] ) );

	// c) ensure valid element width
	if( element_pitch_dbl < element_width_dbl ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidElementPitch", "Element pitch must exceed element width!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 7.) lower bound for temporal frequencies
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 7 ] ) && mxIsScalar( prhs[ 7 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyBound", "Lower frequency bound must be double and scalar!" );
	}

	// b) lower frequency bound
	f_lb_dbl = *( ( double * ) mxGetData( prhs[ 7 ] ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 8.) upper bound for temporal frequencies
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 8 ] ) && mxIsScalar( prhs[ 8 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyBound", "Upper frequency bound must be double and scalar!" );
	}

	// b) upper frequency bound
	f_ub_dbl = *( ( double * ) mxGetData( prhs[ 8 ] ) );

	// c) ensure valid upper frequency bound
	if( f_ub_dbl < f_lb_dbl ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyBound", "Upper frequency bound must equal or exceed the lower frequency bound!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 9.) small-signal sound speed
	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) check mxArray in workspace
	if( !( mxIsDouble( prhs[ 9 ] ) && mxIsScalar( prhs[ 9 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidSoundSpeed", "Average speed of sound must be double and scalar!" );
	}

	// b) small-signal sound speed
	c_0_dbl = *( ( double * ) mxGetData( prhs[ 9 ] ) );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 10.) number of samples added to computed time delay
	//------------------------------------------------------------------------------------------------------------------------------------------
	if( nrhs >= 11 && !mxIsEmpty( prhs[ 10 ] ) && mxIsDouble( prhs[ 10 ] ) )
	{
		N_samples_shift_add_dbl = *( ( double * ) mxGetData( prhs[ 10 ] ) );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 11.) window function for the receive apodization ( object of class windows.window )
	//------------------------------------------------------------------------------------------------------------------------------------------
//TODO: create default F-number
	// a) ensure existence of nonempty window
	// if( nrhs < 12 || mxIsEmpty( prhs[ 11 ] ) )
	// {
	// 	// specify default window (TODO: create MATLAB object)
	// 	index_window_rx = 0;
	// }

	// b) ensure single object of class windows.window
	if( !( isa( prhs[ 11 ], "windows.window" ) && mxIsScalar( prhs[ 11 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidWindow", "Receive window must be a single object of the class windows.window!");
	}

	// c) map classes to window id
	index_window_rx = get_window_id( prhs[ 11 ] );

	// d) check validity of window id
	if( !( ( index_window_rx >= 0 ) && ( index_window_rx < N_WINDOWS ) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidWindow", "Unknown window function!" );

	// e) name of window function
	if( string( prhs[ 11 ], window_rx_name, 128 * sizeof(char) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidWindow", "Error reading window name!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 12.) F-number
	//------------------------------------------------------------------------------------------------------------------------------------------
//TODO: create default F-number
	// a) ensure single object of class f_numbers.f_number
	if( !( isa( prhs[ 12 ], "f_numbers.f_number" ) && mxIsScalar( prhs[ 12 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFNumber", "F-number must be a single object of the class f_numbers.f_number!");
	}

	// b) check if F-number is fixed
	if( isa( prhs[ 12 ], "f_numbers.constant" ) )
	{
		F_number_constant = 1;
	}

	// c) name of F-number
	if( string( prhs[ 12 ], f_number_name, 128 * sizeof(char) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFNumber", "Error reading F-number name!" );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 13.) normalization of the apodization weights
	//------------------------------------------------------------------------------------------------------------------------------------------
	// b) ensure scalar object of class normalizations.normalization
	if( !( isa( prhs[ 13 ], "normalizations.normalization" ) && mxIsScalar( prhs[ 13 ] ) ) )
	{
		mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidNormalization", "Normalization must be a single object of the class normalizations.normalization!");
	}

	// c) map classes to window id
	if( mxIsClass( prhs[ 13 ], "normalizations.windowed" ) )
	{
		index_window_f = get_window_id( mxGetProperty( prhs[ 13 ], 0, "window" ) );
	}

	// d) check validity of window id
	if( !( ( index_window_f >= 0 ) && ( index_window_f < N_WINDOWS ) ) ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:InvalidFrequencyWindow", "Unknown frequency window function!" );

	mexPrintf( "index_window_f = %d\n", index_window_f );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// 14.) index of the GPU
	//------------------------------------------------------------------------------------------------------------------------------------------
	if( nrhs >= 15 && !mxIsEmpty( prhs[ 14 ] ) && mxIsDouble( prhs[ 14 ] ) )
	{
		index_gpu = (int) *( ( double * ) mxGetData( prhs[ 14 ] ) );
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) compute dependent parameters
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) propagation direction of steered plane wave
	//------------------------------------------------------------------------------------------------------------------------------------------
	cos_theta_incident_dbl = cos( theta_incident_dbl );
	sin_theta_incident_dbl = sin( theta_incident_dbl );

	//------------------------------------------------------------------------------------------------------------------------------------------
	// b) array geometry
	//------------------------------------------------------------------------------------------------------------------------------------------
	// compute lateral positions of array elements
	pos_rx_ctr_x = ( float * ) mxMalloc( sizeof( float ) * N_el_rx );
	pos_rx_lb_x = ( float * ) mxMalloc( sizeof( float ) * N_el_rx );
	pos_rx_ub_x = ( float * ) mxMalloc( sizeof( float ) * N_el_rx );

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
	pos_tx_x_ref_dbl = ( cos_theta_incident_dbl >= 0 ) ? - M_el_rx * element_pitch_dbl : M_el_rx * element_pitch_dbl;

	//------------------------------------------------------------------------------------------------------------------------------------------
	// d) compute order of DFT; employ zero-padding to achieve odd order
	//------------------------------------------------------------------------------------------------------------------------------------------
// TODO: determine time duration of focused signal
	// tof_min, tof_max, t_lb, t_ub for all emissions
	if( N_samples_t % 2 == 0 )
	{
		N_order_dft = N_samples_t + 1;
	}
	else
	{
		// N_order_dft = N_samples_t;
		N_order_dft = N_samples_t + 2000;
	}
	// assertion: N_order_dft is odd

	// number of distinct samples in DFT of real-valued signal
	N_samples_dft = N_order_dft / 2 + 1;

	//------------------------------------------------------------------------------------------------------------------------------------------
	// e) compute parameters for bandpass filter
	//------------------------------------------------------------------------------------------------------------------------------------------
	index_f_lb = ceil( f_lb_dbl * N_order_dft / f_s_dbl );
	index_f_ub = floor( f_ub_dbl * N_order_dft / f_s_dbl );

	// check validity of index_f_lb and index_f_ub
	if( ( index_f_lb > index_f_ub ) || ( index_f_lb < 0 ) || ( index_f_ub >= N_samples_dft ) )
	{
		mexPrintf("error: invalid frequency boundaries f_lb and f_ub!");
		return;
	}
	// assertion: index_f_lb <= index_f_ub, index_f_lb >= 0, index_f_ub < N_samples_dft

	// number of non-redundant samples of DFT of order N_order_dft after bandpass filter
	N_samples_dft_bp  = index_f_ub - index_f_lb + 1;

	//------------------------------------------------------------------------------------------------------------------------------------------
	// f) compute frequency-dependent F-number [ in MATLAB ]
	//------------------------------------------------------------------------------------------------------------------------------------------
	// if( F_number_constant )
	// {
	// 	// constant F-number
	// 	f_number_values_dbl = ( double * ) mxGetData( mxGetProperty( prhs[ 12 ], 0, "value" ) );
	// }
	// else
	// {
		// create normalized element pitches in MATLAB workspace
		mxArray* element_pitch_norm_matlab = mxCreateDoubleMatrix( N_samples_dft_bp, 1, mxREAL );

		element_pitch_norm_dbl = (double *) mxGetData( element_pitch_norm_matlab );

		// iterate frequencies
		for( int index_f = 0; index_f < N_samples_dft_bp; index_f++ )
		{

			// normalized element pitch
			element_pitch_norm_dbl[ index_f ] = element_pitch_dbl * ( index_f_lb + index_f ) * f_s_dbl / ( N_order_dft * c_0_dbl );

		}

		// compute values of the frequency-dependent F-number
		f_number_values_dbl = compute_values( prhs[ 12 ], element_pitch_norm_matlab );
	// }

	//------------------------------------------------------------------------------------------------------------------------------------------
	// g) compute frequency-dependent F-number
	//------------------------------------------------------------------------------------------------------------------------------------------
	cufft_batch_size = N_el_rx;

	argument_factor_dbl = 2 * M_PI / N_order_dft;

	// number of blocks to process in parallel
	N_blocks_x = ceil( ( (double) N_pos_lat_x ) / N_THREADS_X );
	N_blocks_z = ceil( ( (double) N_pos_lat_z ) / N_THREADS_Y );

	// determine relevant range of frequency indices
	N_blocks_Omega_bp = ceil( ((double) N_samples_dft_bp) / N_THREADS_PER_BLOCK );

	/*=====================================================================
	 * 2.) print global status information
	 *=====================================================================*/
	mexPrintf( " %s\n", "================================================================================");
	mexPrintf( " GPU BF Toolbox: DAS, steered PW, v.%s (%s), Status Information\n", REVISION, DATE);
	mexPrintf( " (c) 2010-2022, M. F. Schiffner. All rights reserved.\n");
	mexPrintf( " %s\n", "================================================================================");
	mexPrintf( "  %-20s: % 5d x %-5d %5s", "lattice geometry", N_pos_lat_x, N_pos_lat_z, "");
	mexPrintf( " %-22s: %-8d\n", "number of samples", N_samples_t);
	mexPrintf( "  %-20s: %-11d %7s", "number rx channels", N_el_rx, "");
	mexPrintf( " %-22s: %5.1f MHz\n", "sampling rate", f_s_dbl / 1000000);
	mexPrintf( "  %-20s: %7.2f m/s %7s", "average sound speed", c_0_dbl, "");
	mexPrintf( " %-22s: %5.1f MHz (%d)\n", "lower frequency bound", f_lb_dbl / 1000000, index_f_lb);
	mexPrintf( "  %-20s: %6.2f° (%4.1f;%4.1f)", "angle of incidence", theta_incident_dbl * 180 / M_PI, cos_theta_incident_dbl, sin_theta_incident_dbl);
	mexPrintf( " %-22s: %5.1f MHz (%d)\n", "upper frequency bound", f_ub_dbl / 1000000, index_f_ub);
	mexPrintf( "  %-20s: %7.2f mm %8s", "reference position", pos_tx_x_ref_dbl * 1000, "");
	mexPrintf( " %-22s: %-6d\n", "number of frequencies", N_samples_dft_bp);
	mexPrintf( "  %-20s: %-6.2f %12s", "additional shift", N_samples_shift_add_dbl, "");
	mexPrintf( " %-22s: %s (%d)\n", "window", window_rx_name, index_window_rx );
	mexPrintf( "  %-20s: %s\n", "F-number", f_number_name );
	mexPrintf( "  %-20s: %s\n", "normalization", "test" );

    //mexPrintf(" %-22s: %-6.2f um\n", "maximum delta_x", delta_x_max * 1e6);
    //mexPrintf(" %-22s: %-6.2f um\n", "maximum delta_z", delta_z_max * 1e6);

	//=====================================================================
	// 3.) check number of devices supporting CUDA and select device
	//=====================================================================
	int N_devices = 0;         // number of GPUs
// TODO: gpu_bf_set_device
	// check number of devices supporting CUDA
	checkCudaRTErrors( cudaGetDeviceCount( &N_devices ) );
	// mexPrintf( "\nnumber of GPUs supporting CUDA: %d\n", N_devices );
	// mexPrintf("using GPU: %d\n", index_gpu);

	// early exit if no CUDA-enabled GPU was detected
	if( N_devices < 1 ) mexErrMsgIdAndTxt( "gpu_bf_das_pw_rf:NoCUDADevices", "Could not find any CUDA capable devices!" );

	// early exit if selected GPU does not exist
	if( index_gpu >= N_devices || index_gpu < 0 ) mexErrMsgIdAndTxt( "gpu_bf_saft:InvalidCUDADevice", "Invalid GPU selected!" );
	// assertion: N_devices > 0 && 0 <= index_gpu < N_devices

	// set device to operate on
	checkCudaRTErrors( cudaSetDevice( index_gpu ) );

	//=====================================================================
	// 4.) infer and display properties of selected device
	//=====================================================================
	cudaDeviceProp deviceProp;
	int driverVersion = 0;       // driver version
	int runtimeVersion = 0;      // runtime version
	int capability_double = 0;   // flag for capability for double precision

	// check properties of chosen GPU
	checkCudaRTErrors( cudaGetDeviceProperties( &deviceProp, index_gpu ) );

	// check capabilities for double precision
	if( deviceProp.major == 1 && deviceProp.minor >= 3 )
	{
		capability_double = 1;
	}
	if( deviceProp.major >= 2)
	{
		capability_double = 1;
	}

	// get driver version
	checkCudaRTErrors( cudaDriverGetVersion( &driverVersion ) );

	// get runtime version
	checkCudaRTErrors( cudaRuntimeGetVersion( &runtimeVersion ) );

	/*=====================================================================
     * 5.) print GPU status information
     *=====================================================================*/
    mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
    mexPrintf(" Information for GPU device %-1d of %-1d:\n", index_gpu, N_devices);
    mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
    mexPrintf("  %-20s: %-19s", "name of device", deviceProp.name, "");
    mexPrintf(" %-22s: %-3d\n", "num. multiprocessors", deviceProp.multiProcessorCount);
    mexPrintf("  %-20s: %-1d.%-1d %15s", "compute capability", deviceProp.major, deviceProp.minor, "");
    mexPrintf(" %-22s: %-8d\n", "double precision", capability_double);
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

	/*=====================================================================
	 * 6.) start calculations depending on selected precision
	 *=====================================================================*/
	if(0)
	{
		mexPrintf("\nprecision of calculations: double\n");
	}
	else
	{

		/*===============================================================
		 * compute required size for arrays in bytes (FLOAT)
		 *===============================================================*/
		size_bytes_positions_x = sizeof( t_float_gpu ) * N_pos_lat_x;
		size_bytes_positions_z = sizeof( t_float_gpu ) * N_pos_lat_z;
		size_bytes_RF_data = sizeof( cufftComplex ) * N_samples_dft * N_el_rx;
		size_bytes_signal  = sizeof( t_float_gpu ) * N_samples_t;
		size_bytes_positions_rx = sizeof( t_float_gpu ) * N_el_rx;
		size_bytes_weights = sizeof( t_float_gpu ) * N_pos_lat_x * N_pos_lat_z;
		size_bytes_f_numbers = sizeof( float ) * N_samples_dft_bp;

		size_bytes_DAS_image_RF_complex = sizeof( cufftComplex ) * N_pos_lat_x * N_pos_lat_z;

		//===============================================================
		// allocate memory for the MATLAB output
		//===============================================================
		// complex-valued DAS image
		plhs[ 0 ] = mxCreateNumericMatrix( N_pos_lat_z, N_pos_lat_x, mxSINGLE_CLASS, mxCOMPLEX );
		image_DAS_real = (float*) mxGetData( plhs[ 0 ] );
		image_DAS_imag = (float*) mxGetImagData( plhs[ 0 ] );

		// total weights for grid points (host)
		plhs[ 1 ] = mxCreateNumericMatrix( N_pos_lat_z, N_pos_lat_x, mxSINGLE_CLASS, mxREAL );
		h_weights = (float*) mxGetData( plhs[ 1 ] );

		/*===============================================================
		 * convert input data to float
		 *===============================================================*/
		pos_lat_x = (float *) mxMalloc( size_bytes_positions_x );
		pos_lat_z = (float *) mxMalloc( size_bytes_positions_z );
		data_RF   = (float *) mxMalloc( sizeof( float ) * N_samples_t * N_el_rx );
		f_number_values_flt = (float *) mxMalloc( size_bytes_f_numbers );

		// convert lateral image coordinates
		for( l_x = 0; l_x < N_pos_lat_x; l_x++ )
		{
   			pos_lat_x[ l_x ] = (float) pos_lat_x_dbl[ l_x ];
		}

		// convert axial image coordinates
		for( l_z = 0; l_z < N_pos_lat_z; l_z++ )
		{
			pos_lat_z[ l_z ] = (float) pos_lat_z_dbl[ l_z ];
		}

		// convert RF data
		for( k_rx = 0; k_rx < N_el_rx; k_rx++ )
		{
			for( l_z = 0; l_z < N_samples_t; l_z++ )
			{
				data_RF[ k_rx * N_samples_t + l_z ] = (float) data_RF_dbl[ k_rx * N_samples_t + l_z ];      
			}
		}

		// convert frequency-dependent F-number
		for( int index_f = 0; index_f < N_samples_dft_bp; index_f++ )
		{
			f_number_values_flt[ index_f ] = (float) f_number_values_dbl[ index_f ];
			// mexPrintf( "f_number_values_flt[ %d ] = %.2f\n", index_f, f_number_values_flt[ index_f ] );
		}
		// mxDestroyArray( lhs[ 0 ] );

		// convert components of direction vector
		cos_theta_incident_flt = (float) cos_theta_incident_dbl;
		sin_theta_incident_flt = (float) sin_theta_incident_dbl;

		// convert reference position
		pos_tx_ctr_x_ref_flt = (float) pos_tx_x_ref_dbl;

		// convert additional time shift
		N_samples_shift_add_flt = (float) N_samples_shift_add_dbl;

		// convert sampling rate and speed of sound
		f_s_over_c_0_flt = (float) (f_s_dbl / c_0_dbl);

		// convert factor in exponent
		argument_factor_flt = (float) argument_factor_dbl;

		/*===============================================================
		 * print memory status information
		 *===============================================================*/
		size_bytes_total = size_bytes_positions_x + size_bytes_positions_z + size_bytes_RF_data + 3 * size_bytes_positions_rx + size_bytes_f_numbers + size_bytes_DAS_image_RF_complex + size_bytes_weights;

		// get available amount of global memory on GPU device
		checkCudaRTErrors( cudaMemGetInfo( &mem_available, &(deviceProp.totalGlobalMem) ) );

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

		/*===============================================================
		 * allocate memory on device
		 *===============================================================*/
		// a) grid points
		checkCudaRTErrors( cudaMalloc( (void **) &d_pos_lat_x_flt, size_bytes_positions_x ) );
		checkCudaRTErrors( cudaMalloc( (void **) &d_pos_lat_z_flt, size_bytes_positions_z ) );

		// b) allocate memory for RF data and its DFT
		checkCudaRTErrors( cudaMalloc( (void **) &data_RF_flt_gpu, size_bytes_RF_data ) );

		// initialize RF data (necessary because of zero-padding)
		checkCudaRTErrors( cudaMemset( data_RF_flt_gpu, 0, size_bytes_RF_data ) );

		// c) array geometry
		checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_ctr_x, size_bytes_positions_rx) );
		checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_lb_x, size_bytes_positions_rx) );
		checkCudaRTErrors( cudaMalloc( (void **) &d_pos_rx_ub_x, size_bytes_positions_rx) );

		// d) allocate memory for frequency-dependent F-number
		checkCudaRTErrors( cudaMalloc( (void **) &d_f_number_values, size_bytes_f_numbers ) );

		// e) allocate memory for DAS image
		checkCudaRTErrors( cudaMalloc( (void **) &d_image_flt_cpx, size_bytes_DAS_image_RF_complex ) );

		// initialize complex-valued DAS image
		checkCudaRTErrors( cudaMemset( d_image_flt_cpx, 0, size_bytes_DAS_image_RF_complex ) );

		// f) allocate memory for weights
		checkCudaRTErrors( cudaMalloc( (void**) &d_weights, size_bytes_weights ) );

		// initialize weights
		checkCudaRTErrors( cudaMemset( d_weights, 0, size_bytes_weights ) );

		/*=================================================================
		 * 1.) copy data to device
		 *=================================================================*/
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf(" running algorithm:\n");
		mexPrintf(" %s\n", "--------------------------------------------------------------------------------");
		mexPrintf(" %4s 1.) %s... ", "", "transfering data from host to device");

		// create events for time-measurement
		checkCudaRTErrors( cudaEventCreate( &start ) );
		checkCudaRTErrors( cudaEventCreate( &stop ) );

		// place start event into the default stream
		checkCudaRTErrors( cudaEventRecord( start, 0 ) );

		// a) grid points
		checkCudaRTErrors( cudaMemcpy( d_pos_lat_x_flt, pos_lat_x, size_bytes_positions_x, cudaMemcpyHostToDevice ) );
		checkCudaRTErrors( cudaMemcpy( d_pos_lat_z_flt, pos_lat_z, size_bytes_positions_z, cudaMemcpyHostToDevice ) );

		// b) correctly pad RF data to allow in place transform
		for( k_rx = 0; k_rx < N_el_rx; k_rx++ )
		{
			checkCudaRTErrors( cudaMemcpy( data_RF_flt_gpu + k_rx * N_samples_dft * 2, &( data_RF[ k_rx * N_samples_t ] ), size_bytes_signal, cudaMemcpyHostToDevice ) );
		}

		// c) array geometry
		checkCudaRTErrors( cudaMemcpy( d_pos_rx_ctr_x, pos_rx_ctr_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );
		checkCudaRTErrors( cudaMemcpy( d_pos_rx_lb_x, pos_rx_lb_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );
		checkCudaRTErrors( cudaMemcpy( d_pos_rx_ub_x, pos_rx_ub_x, size_bytes_positions_rx, cudaMemcpyHostToDevice ) );

		// d) frequency-dependent F-number
		checkCudaRTErrors( cudaMemcpy( d_f_number_values, f_number_values_flt, size_bytes_f_numbers, cudaMemcpyHostToDevice ) );

		// place stop event into the default stream
		checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
		checkCudaRTErrors( cudaEventSynchronize( stop ) );

		// compute elapsed time
		checkCudaRTErrors( cudaEventElapsedTime( &time_transfer_to_device, start, stop ) );

		mexPrintf("done! (%.3f ms)\n", time_transfer_to_device);

		// flush cache
		mexEvalString( "drawnow;" );

		//=====================================================================
		// 2.) compute DFTs of data_RF along first MATLAB dimension; choose appropriate order; in place
		//=====================================================================
		mexPrintf(" %4s 2.) %s... ", "", "computing DFT");

		// place start event into the default stream
		checkCudaRTErrors( cudaEventRecord(start, 0) );

		cufftHandle plan;

		// create plan for DFT
		checkCudaFFTErrors( cufftPlan1d( &plan, N_order_dft, CUFFT_R2C, cufft_batch_size ) );

		// execute DFT plan
		checkCudaFFTErrors( cufftExecR2C( plan, (cufftReal *) data_RF_flt_gpu, (cufftComplex *) data_RF_flt_gpu ) );

		// destroy plan after execution
		checkCudaFFTErrors( cufftDestroy( plan ) );

		// place stop event into the default stream
		checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
		checkCudaRTErrors( cudaEventSynchronize( stop ) );

		// compute elapsed time
		checkCudaRTErrors( cudaEventElapsedTime( &time_dft, start, stop ) );

		mexPrintf("done! (%.3f ms)\n", time_dft);

		// flush cache
		mexEvalString( "drawnow;" );

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// 3.) copy addresses of __device__ window functions to host
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		checkCudaRTErrors( cudaMemcpyFromSymbol( &h_windows, d_windows, N_WINDOWS * sizeof( t_window_ptr ) ) );

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// 4.) compute DAS image on GPU
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		mexPrintf(" %4s 3.) %s... ", "", "computing DAS image on GPU");

		dim3 threadsPerBlock(N_THREADS_X, N_THREADS_Y);  /* use max number of threads per block for additions*/
		dim3 numBlocks(N_blocks_x, N_blocks_z);

		// place start event into the default stream
		checkCudaRTErrors( cudaEventRecord( start, 0 ) );

		// iterate rx channels
		for( k_rx = 0; k_rx < N_el_rx; k_rx++ )
		{

			// compute DAS image for current tx/rx configuration
			if( F_number_constant )
			{
				das_F_number_constant<<<numBlocks, threadsPerBlock>>>( d_image_flt_cpx, d_weights,
																	   ((cufftComplex *) data_RF_flt_gpu) + k_rx * N_samples_dft,
																	   d_pos_lat_x_flt, d_pos_lat_z_flt,
																	   d_pos_rx_ctr_x, d_pos_rx_lb_x, d_pos_rx_ub_x,
																	   ( t_float_gpu ) f_number_values_dbl[ 0 ],
																	   h_windows[ index_window_rx ],
																	   argument_factor_flt, N_blocks_Omega_bp, index_f_lb, index_f_ub,
																	   N_pos_lat_x, N_pos_lat_z, cos_theta_incident_flt, sin_theta_incident_flt, pos_tx_ctr_x_ref_flt,
																	   N_samples_shift_add_flt, f_s_over_c_0_flt, k_rx,
																	   ( t_float_gpu ) element_pitch_dbl, ( t_float_gpu ) M_el_rx, N_el_rx );
			}
			else
			{
				if( index_window_f == 0)
				{
					// boxcar window
					das_F_number_frequency_dependent<<<numBlocks, threadsPerBlock>>>( d_image_flt_cpx,
																					  ((cufftComplex *) data_RF_flt_gpu) + k_rx * N_samples_dft,
																					  d_pos_lat_x_flt, d_pos_lat_z_flt,
																					  d_pos_rx_ctr_x, d_pos_rx_lb_x, d_pos_rx_ub_x,
																					  d_f_number_values,
																					  h_windows[ index_window_rx ],
																					  argument_factor_flt, N_blocks_Omega_bp, index_f_lb, index_f_ub,
																					  N_pos_lat_x, N_pos_lat_z, cos_theta_incident_flt, sin_theta_incident_flt, pos_tx_ctr_x_ref_flt,
																					  N_samples_shift_add_flt, f_s_over_c_0_flt, k_rx,
																					  ( t_float_gpu ) element_pitch_dbl, ( t_float_gpu ) M_el_rx, N_el_rx );
				}
				else
				{
					// other window
					das_F_number_frequency_dependent_window<<<numBlocks, threadsPerBlock>>>( d_image_flt_cpx,
																							 ((cufftComplex *) data_RF_flt_gpu) + k_rx * N_samples_dft,
																							 d_pos_lat_x_flt, d_pos_lat_z_flt,
																							 d_pos_rx_ctr_x, d_pos_rx_lb_x, d_pos_rx_ub_x,
																							 d_f_number_values,
																							 h_windows[ index_window_rx ], h_windows[ index_window_f ],
																							 argument_factor_flt, N_blocks_Omega_bp, index_f_lb, N_samples_dft_bp,
																							 N_pos_lat_x, N_pos_lat_z, cos_theta_incident_flt, sin_theta_incident_flt, pos_tx_ctr_x_ref_flt,
																							 N_samples_shift_add_flt, f_s_over_c_0_flt, k_rx,
																							 ( t_float_gpu ) element_pitch_dbl, ( t_float_gpu ) M_el_rx, N_el_rx );
				}
			} // if( F_number_constant )

		} // for( k_rx = 0; k_rx < N_el_rx; k_rx++ )

		// place stop event into the default stream
		checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
		checkCudaRTErrors( cudaEventSynchronize( stop ) );

		// compute elapsed time
		checkCudaRTErrors( cudaEventElapsedTime( &time_kernel, start, stop ) );

		mexPrintf("done! (%.3f ms)\n", time_kernel);

		// flush cache
		mexEvalString( "drawnow;" );

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// 5.) copy DAS image and weights to host
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// allocate memory for results
		image_DAS_flt_cpx = (cufftComplex *) mxMalloc( size_bytes_DAS_image_RF_complex );

		mexPrintf(" %4s 4.) %s... ", "", "transfering data to host");

		// place start event into the default stream
		checkCudaRTErrors( cudaEventRecord( start, 0 ) );

		// copy results to host
		checkCudaRTErrors( cudaMemcpy( image_DAS_flt_cpx, d_image_flt_cpx, size_bytes_DAS_image_RF_complex, cudaMemcpyDeviceToHost ) );

		// copy weights to host
		if( F_number_constant )
		{
			checkCudaRTErrors( cudaMemcpy( h_weights, d_weights, size_bytes_weights, cudaMemcpyDeviceToHost ) );
		}

		// extract complex-valued pixels and divide by weights (if appropriate)
		for( l_z = 0; l_z < N_pos_lat_z; l_z++ )
		{
			for( l_x = 0; l_x < N_pos_lat_x; l_x++ )
			{

				image_DAS_real[ l_x * N_pos_lat_z + l_z ] = image_DAS_flt_cpx[ l_x * N_pos_lat_z + l_z ].x;
				image_DAS_imag[ l_x * N_pos_lat_z + l_z ] = image_DAS_flt_cpx[ l_x * N_pos_lat_z + l_z ].y;

				// normalize pixels for constant F-number
				if( F_number_constant )
				{
					// check how many rx channels contributed to current pixel
					if( h_weights[ l_x * N_pos_lat_z + l_z ] > 0 )
					{
						image_DAS_real[ l_x * N_pos_lat_z + l_z ] = image_DAS_real[ l_x * N_pos_lat_z + l_z ] / ( N_order_dft * h_weights[ l_x * N_pos_lat_z + l_z ] );
						image_DAS_imag[ l_x * N_pos_lat_z + l_z ] = image_DAS_imag[ l_x * N_pos_lat_z + l_z ] / ( N_order_dft * h_weights[ l_x * N_pos_lat_z + l_z ] );
					}
				} // if( F_number_constant )

				// else
				// {
				// 	// no rx channel contributed to current pixel
				// 	image_DAS_real[ l_x * N_pos_lat_z + l_z ] = 0;
				// 	image_DAS_imag[ l_x * N_pos_lat_z + l_z ] = 0;
				// }
			} // for( l_x = 0; l_x < N_pos_lat_x; l_x++ )
		} // for( l_z = 0; l_z < N_pos_lat_z; l_z++ )

		// place stop event into the default stream
		checkCudaRTErrors( cudaEventRecord( stop, 0 ) );
		checkCudaRTErrors( cudaEventSynchronize( stop ) );

		// compute elapsed time
		checkCudaRTErrors( cudaEventElapsedTime( &time_transfer_to_host, start, stop ) );

		mexPrintf("done! (%.3f ms)\n", time_transfer_to_host);

		// total elapsed time
		time_total = time_transfer_to_device + time_dft + time_kernel + time_transfer_to_host;

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// 6.) print performance statistics
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// 7.) clean up memory
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//------------------------------------------------------------------------------------------------------------------------------------------
		// a) device memory
		//------------------------------------------------------------------------------------------------------------------------------------------
		checkCudaRTErrors( cudaFree( d_f_number_values ) );
		checkCudaRTErrors( cudaFree( data_RF_flt_gpu ) );
		checkCudaRTErrors( cudaFree( d_image_flt_cpx ) );
		checkCudaRTErrors( cudaFree( d_weights ) );
		checkCudaRTErrors( cudaFree( d_pos_lat_x_flt ) );
		checkCudaRTErrors( cudaFree( d_pos_lat_z_flt ) );
		checkCudaRTErrors( cudaFree( d_pos_rx_ctr_x ) );
		checkCudaRTErrors( cudaFree( d_pos_rx_lb_x ) );
		checkCudaRTErrors( cudaFree( d_pos_rx_ub_x ) );
		checkCudaRTErrors( cudaEventDestroy( start ) );
		checkCudaRTErrors( cudaEventDestroy( stop ) );

		//------------------------------------------------------------------------------------------------------------------------------------------
		// b) host memory
		//------------------------------------------------------------------------------------------------------------------------------------------
		mxFree( pos_lat_x );
		mxFree( pos_lat_z );
		mxFree( data_RF );
		mxFree( f_number_values_flt );
		mxFree( pos_rx_ctr_x );
		mxFree( pos_rx_lb_x  );
		mxFree( pos_rx_ub_x  );
		mxFree( image_DAS_flt_cpx );
	}

} // void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
