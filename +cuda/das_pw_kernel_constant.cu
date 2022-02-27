//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// author: Martin F. Schiffner
// date: 2022-02-20
// modified: 2022-02-20

// __device__ t_float_gpu distance_pw( t_float_gpu pos_focus_x, t_float_gpu pos_focus_z, t_float_gpu pos_el_ctr_x )
// {

// }

__global__ void das_F_number_constant( t_float_complex_gpu* image, t_float_gpu* weights,
									   const t_float_complex_gpu* data_RF_dft, const t_float_gpu* pos_lat_x, const t_float_gpu* pos_lat_z,
									   const t_float_gpu* pos_rx_ctr_x, const t_float_gpu* pos_rx_lb_x, const t_float_gpu* pos_rx_ub_x,
									   const t_float_gpu f_number_value, t_window_ptr apo_function_rx,
									   t_float_gpu argument_factor_flt, int N_blocks_Omega_bp, int index_f_lb, int index_f_ub,
									   int N_pos_lat_x, int N_pos_lat_z, t_float_gpu cos_theta_incident, t_float_gpu sin_theta_incident, t_float_gpu pos_tx_ctr_x_ref,
									   t_float_gpu N_samples_shift_add, t_float_gpu f_s_over_c_0, int index_rx,
									   t_float_gpu element_pitch, t_float_gpu M_elements, int N_elements )
{

	// each thread block computes one sub-matrix of final DAS image

	// initialize local variables
	int l_x = blockIdx.x * blockDim.x + threadIdx.x;  //current x-index
	int l_z = blockIdx.y * blockDim.y + threadIdx.y;  //current z-index
	int index_thread = threadIdx.y * blockDim.x + threadIdx.x;

	int index_f = 0; int index_f_lb_block = 0; int index_f_in_block = 0;

	cufftComplex pixel_value;
	pixel_value.x = 0;
	pixel_value.y = 0;

	// compute time-shift required for current x and z coordinates
	// compute argument of twiddle factor
	float pos_focus_x = 0.0f;
	float distance_x = 0.0f; float distance_z = 0.0f;
	float distance_sw = 0.0f;
	float argument = 0.0f;

	float width_aperture_desired_over_two = 0.0f;
	float width_aperture_over_two_max = element_pitch * N_elements;
	float N_samples_shift = 0.0f;
	float apodization = 0.0f;

	float argument_act = 0.0f;
	float value_cosine = 0.0f; float value_sine = 0.0f;
	int index_aperture = 0;

	// shared memory
	__shared__ cufftComplex data_RF_dft_times_two_shared[ N_THREADS_PER_BLOCK ];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// compute round-trip time-of-flight and argument of complex exponential function
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if( l_x < N_pos_lat_x && l_z < N_pos_lat_z )
	{

		// compute round-trip time-of-flight
		pos_focus_x = pos_lat_x[ l_x ];
		distance_x = pos_rx_ctr_x[ index_rx ] - pos_focus_x;
		distance_z = pos_lat_z[ l_z ];
		distance_sw = __fsqrt_rn( distance_x * distance_x + distance_z * distance_z );
		distance_sw = distance_sw + cos_theta_incident * ( pos_focus_x - pos_tx_ctr_x_ref ) + sin_theta_incident * distance_z;

		// argument for complex exponential in inverse DFT
		//N_samples_shift = roundf( f_s_over_c_0 * distance_sw  ) + N_samples_shift_add;
		N_samples_shift = f_s_over_c_0 * distance_sw + N_samples_shift_add;
		argument = argument_factor_flt * N_samples_shift;

		// compute apodization weight for current frequency using F-number
		if( f_number_value > FLT_EPSILON )
		{
			// desired half-width of the receive aperture
			width_aperture_desired_over_two = distance_z / ( 2 * f_number_value );
		}
		else // f_number_value <= FLT_EPSILON
		{
			// full aperture
			width_aperture_desired_over_two = width_aperture_over_two_max;
		}

		// check location of array element
		if( distance_x < 0 )
		{
			// array element belongs to left aperture
			index_aperture = ( int ) fmaxf( ceilf( M_elements + ( pos_focus_x - width_aperture_desired_over_two ) / element_pitch ), 0 );
			apodization = apply( apo_function_rx, distance_x, pos_focus_x - pos_rx_lb_x[ index_aperture ] );
		}
		else
		{
			// array element belongs to right aperture
			index_aperture = ( int ) fminf( floorf( M_elements + ( pos_focus_x + width_aperture_desired_over_two ) / element_pitch ), N_elements - 1 );
			apodization = apply( apo_function_rx, distance_x, pos_rx_ub_x[ index_aperture ] - pos_focus_x );
		}

	} // if( l_x < N_pos_lat_x && l_z < N_pos_lat_z )

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// loop over all frequency blocks that are required to compute sub-matrix of final image
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for( int index_block_f = 0; index_block_f < N_blocks_Omega_bp; index_block_f++ )
	{

		// each thread in current block transfers one complex DFT value from device memory to shared memory

		// compute current frequency index, take into account thread ID
		index_f_lb_block = index_f_lb + index_block_f * N_THREADS_PER_BLOCK;
		index_f = index_f_lb_block + index_thread;

		// check validity of frequency
		if( index_f <= index_f_ub )
		{
			// transfer complex DFT value from device memory to shared memory
			// each thread loads the value for one frequency index
			data_RF_dft_times_two_shared[ index_thread ].x = 2 * data_RF_dft[ index_f ].x;
			data_RF_dft_times_two_shared[ index_thread ].y = 2 * data_RF_dft[ index_f ].y;
		}
		else
		{
			data_RF_dft_times_two_shared[ index_thread ].x = 0.0f;
			data_RF_dft_times_two_shared[ index_thread ].y = 0.0f;
		}

		// synchronize threads to ensure a complete data set before starting calculations
		__syncthreads();

		// special case: zero frequency
		if( ( index_f_lb == 0 ) && ( index_block_f == 0 ) )
		{
			pixel_value.x = - 0.5f * data_RF_dft_times_two_shared[ 0 ].x;
			pixel_value.y = - data_RF_dft_times_two_shared[ 0 ].y;
		}

		// iterate frequencies in current block
		for( index_f_in_block = 0; index_f_in_block < N_THREADS_PER_BLOCK; index_f_in_block++ )
		{

			// argument of complex exponential
			argument_act = argument * ( index_f_lb_block + index_f_in_block );
			__sincosf( argument_act, &value_sine, &value_cosine );

			// assign pixel values
			pixel_value.x = pixel_value.x + data_RF_dft_times_two_shared[ index_f_in_block ].x * value_cosine - data_RF_dft_times_two_shared[ index_f_in_block ].y * value_sine;
			pixel_value.y = pixel_value.y + data_RF_dft_times_two_shared[ index_f_in_block ].x * value_sine + data_RF_dft_times_two_shared[ index_f_in_block ].y * value_cosine;

		} // for( index_f_in_block = 0; index_f_in_block < N_THREADS_PER_BLOCK; index_f_in_block++ )

		// synchronize threads to ensure complete calculations before loading new frequency data
		__syncthreads();

	} // for( index_block_f = 0; index_block_f < N_blocks_Omega_bp; index_block_f++ )

	if( l_x < N_pos_lat_x && l_z < N_pos_lat_z && apodization > FLT_EPSILON )
	{
		image[ l_x * N_pos_lat_z + l_z ].x += apodization * pixel_value.x;
		image[ l_x * N_pos_lat_z + l_z ].y += apodization * pixel_value.y;

		weights[ l_x * N_pos_lat_z + l_z ] += apodization;
	}

}