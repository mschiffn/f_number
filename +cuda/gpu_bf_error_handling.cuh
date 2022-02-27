//-------------------------------------------------------------------------
// canonical error handling (macros and function prototypes)
//-------------------------------------------------------------------------
#ifndef __GPU_BF_ERROR_HANDLING_CUH__
#define __GPU_BF_ERROR_HANDLING_CUH__

#include <cuda.h>
#include <cufft.h>

// wrapper macros
#define checkCudaRTErrors( val ) _checkCudaRTErrors( (val), #val, __FILE__, __LINE__ )
#define checkCudaFFTErrors( val ) _checkCudaFFTErrors( (val), #val, __FILE__, __LINE__ )

// error handling
void _checkCudaRTErrors( const cudaError_t result, const char* const str_command, const char* const str_filename, const int index_line );
void _checkCudaFFTErrors( const cufftResult_t result, const char* const str_command, const char* const str_filename, const int index_line );

#endif
