// load required header files
#include <mex.h>
#include <cuda.h>
#include <cufft.h>

// define data types
typedef float t_float_gpu;
typedef cufftComplex t_float_complex_gpu;

// window function pointer type
// typedef t_float_gpu (*t_window_ptr)( const t_float_gpu, const t_float_gpu );

/*=======================================================================
 * define function prototypes
 *=======================================================================*/

// device kernels
// __global__ void DAS_block_kernel( t_float_complex_gpu* image, t_float_gpu* weights,
//                                   const t_float_complex_gpu* data_RF_dft, const t_float_gpu* pos_lat_x, const t_float_gpu* pos_lat_z,
//                                   const t_float_gpu* pos_rx_ctr_x, const t_float_gpu* pos_rx_lb_x, const t_float_gpu* pos_rx_ub_x,
//                                   const t_float_gpu* f_number_values, t_window_ptr apo_function_rx,
//                                   t_float_gpu argument_factor, int N_blocks_Omega_bp, int index_Omega_lb, int index_Omega_ub,
//                                   int N_pos_lat_x, int N_pos_lat_z, t_float_gpu cos_theta_incident, t_float_gpu sin_theta_incident, t_float_gpu pos_tx_ctr_x_ref,
//                                   t_float_gpu N_samples_shift_add, int N_order_dft, t_float_gpu f_s_over_c_0, int index_rx,
//                                   t_float_gpu element_pitch, t_float_gpu M_elements, int N_elements );

// device functions
// __device__ t_float_gpu boxcar( const t_float_gpu l, const t_float_gpu L_over_2 );
// __device__ t_float_gpu hanning( const t_float_gpu l, const t_float_gpu L_over_2 );
