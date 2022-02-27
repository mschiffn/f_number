//-------------------------------------------------------------------------
// canonical error handling (implementation)
//-------------------------------------------------------------------------
#include "gpu_bf_error_handling.cuh"

//-------------------------------------------------------------------------
// CUDA runtime API
//-------------------------------------------------------------------------
void _checkCudaRTErrors( const cudaError_t result, const char* const str_command, const char* const str_filename, int const index_line )
{
	if( cudaSuccess != result )
	{

		// print error location and CUDA error strings
		mexPrintf( "\nCUDA error in line %d of file \"%s\":\n", index_line, str_filename );
		mexPrintf( "\t\"%s\"\n", str_command );
		mexPrintf( "\tcode = %d (%s): %s\n", static_cast<unsigned int>( result ), cudaGetErrorName( result ), cudaGetErrorString( result ) );

		// reset device to clean memory before exit
		cudaDeviceReset();

		// print error message, exit program
		mexErrMsgIdAndTxt( "FMM:ErrorCUDA", "CUDA error!" );

	}
} // void _checkCudaRTErrors( const cudaError_t result, const char* const str_command, const char* const str_filename, int const index_line )

//-------------------------------------------------------------------------
// cuFFT library
//-------------------------------------------------------------------------
void _checkCudaFFTErrors( const cufftResult_t result, const char* const str_command, const char* const str_filename, int const index_line )
{
	if( CUFFT_SUCCESS != result )
	{

		// print error location and CUDA error strings
		mexPrintf( "CUDA error at %s:%d:\n", str_filename, index_line );
		mexPrintf( "\t\"%s\"\n", str_command );
		mexPrintf( "code = %d\n", static_cast<unsigned int>( result ) );

		// reset device to clean memory before exit
		cudaDeviceReset();

		// print error message, exit program
		mexErrMsgIdAndTxt( "FMM:ErrorCUDA", "CUDA error!" );

	}
} // void _checkCudaFFTErrors( const cudaError_t result, const char* const str_command, const char* const str_filename, int const index_line )
