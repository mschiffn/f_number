//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// DAS kernel, steered PW, frequency-dependent F-number
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// author: Martin F. Schiffner
// date: 2022-02-20
// modified: 2023-12-14

__global__ void das_F_number_frequency_dependent( t_float_complex_gpu* image,
													const t_float_complex_gpu* data_RF_dft,
													const t_float_gpu* positions_x, const t_float_gpu* positions_z,
													const t_float_gpu* pos_rx_ctr_x, const t_float_gpu* pos_rx_lb_x, const t_float_gpu* pos_rx_ub_x,
													const t_float_gpu* f_number_values, t_window_ptr window_rx,
													const t_float_gpu argument_factor, int N_blocks_Omega_bp, int index_f_lb, int index_f_ub,
													int N_pos_lat_x, int N_pos_lat_z, t_float_gpu e_steering_x, t_float_gpu e_steering_z, t_float_gpu pos_tx_ctr_x_ref,
													t_float_gpu argument_add, const int index_rx,
													t_float_gpu element_pitch, t_float_gpu M_elements, int N_elements )
{

	// each thread block computes one sub-matrix of final DAS image
// TODO: load element_width and replace pos_rx_lb_x and pos_rx_ub_x

	// initialize local variables
	int l_x = blockIdx.x * blockDim.x + threadIdx.x;  //current x-index
	int l_z = blockIdx.y * blockDim.y + threadIdx.y;  //current z-index
	int index_thread = threadIdx.y * blockDim.x + threadIdx.x;

	int index_f_shift = 0;
	int index_f = 0; int index_f_in_band = 0; int index_f_in_block = 0;

	cufftComplex pixel_value;
	pixel_value.x = 0;
	pixel_value.y = 0;

	// compute time-shift required for current x and z coordinates
	// compute argument of twiddle factor
	float pos_focus_x = 0.0f;
	float distance_x = 0.0f; float distance_z = 0.0f;
	float argument = 0.0f;

	float width_aperture_left_over_two = 0.0f; float width_aperture_right_over_two = 0.0f;
	float width_aperture_over_two_max = element_pitch * N_elements;
	float apodization = 1.0f; float apodization_act = 1.0f;	float apodization_sum = 0.0f;

	float distance_x_act = 0.0f;
	float argument_act = 0.0f;
	float value_cosine = 0.0f; float value_sine = 0.0f;
	int index_aperture_lb = 0; int index_aperture_ub = 0;

	// shared memory
	__shared__ float f_number_values_shared[ N_THREADS_PER_BLOCK ];
	__shared__ cufftComplex data_RF_dft_times_two_shared[ N_THREADS_PER_BLOCK ];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) compute argument of complex exponential function
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if( l_x < N_pos_lat_x && l_z < N_pos_lat_z )
	{

		//--------------------------------------------------------------------------------------------------------------------------------------
		// 1.) compute argument
		//--------------------------------------------------------------------------------------------------------------------------------------
		// a) compute round-trip distance
		pos_focus_x = positions_x[ l_x ];
		distance_x = pos_focus_x - pos_rx_ctr_x[ index_rx ];
		distance_z = positions_z[ l_z ];
		float distance_sw = e_steering_x * ( pos_focus_x - pos_tx_ctr_x_ref ) + e_steering_z * distance_z;
		distance_sw = distance_sw + __fsqrt_rn( distance_x * distance_x + distance_z * distance_z );

		// b) argument for complex exponential in inverse DFT
		argument = argument_factor * distance_sw + argument_add;

	} // if( l_x < N_pos_lat_x && l_z < N_pos_lat_z )

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 2.) iterate all frequency blocks that are required to compute sub-matrix of final image
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for( int index_block_f = 0; index_block_f < N_blocks_Omega_bp; index_block_f++ )
	{

		// each thread in current block transfers one complex DFT value from device memory to shared memory

		// compute current frequency index, take into account thread ID
		index_f_shift = index_block_f * N_THREADS_PER_BLOCK;
		index_f_in_band = index_f_shift + index_thread;
		index_f = index_f_lb + index_f_in_band;

		// check validity of frequency
		if( index_f <= index_f_ub )
		{
			// transfer complex DFT value from device memory to shared memory
			// each thread loads the value for one frequency index
			data_RF_dft_times_two_shared[ index_thread ].x = 2 * data_RF_dft[ index_f ].x;
			data_RF_dft_times_two_shared[ index_thread ].y = 2 * data_RF_dft[ index_f ].y;

			f_number_values_shared[ index_thread ] = f_number_values[ index_f_in_band ];
		}
		else
		{
			data_RF_dft_times_two_shared[ index_thread ].x = 0.0f;
			data_RF_dft_times_two_shared[ index_thread ].y = 0.0f;

			f_number_values_shared[ index_thread ] = 0.0f;
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

			// current frequency index relative to index_f_lb
			index_f_in_band = index_f_shift + index_f_in_block;

			// check voxel validity
			if( l_x < N_pos_lat_x && l_z < N_pos_lat_z )
			{
				// compute apodization weight for current frequency using F-number
				// boundaries of the receive aperture
				// a) default bounds (F-number of zero)
				index_aperture_lb = 0;
				index_aperture_ub = N_elements - 1;

				// b) bounds for positive F-number
				if( f_number_values_shared[ index_f_in_block ] > FLT_EPSILON )
				{
					// desired half-width of the receive aperture
					float width_aperture_over_two_desired = distance_z / ( 2 * f_number_values_shared[ index_f_in_block ] );
					index_aperture_lb = ( int ) fmaxf( ceilf( M_elements + ( pos_focus_x - width_aperture_over_two_desired ) / element_pitch ), 0 );
					index_aperture_ub = ( int ) fminf( floorf( M_elements + ( pos_focus_x + width_aperture_over_two_desired ) / element_pitch ), N_elements - 1 );
				}

				// next frequency if rx element is outside the receive aperture
				if( ( index_rx < index_aperture_lb ) || ( index_rx > index_aperture_ub ) ) continue;
				// assertion: index_rx >= index_aperture_lb && index_rx <= index_aperture_ub

				// half-widths of the receive aperture
				width_aperture_left_over_two = pos_focus_x - pos_rx_lb_x[ index_aperture_lb ];
				width_aperture_right_over_two = pos_rx_ub_x[ index_aperture_ub ] - pos_focus_x;

				// reset sum of apodization weights
				apodization_sum = 0.0f;

				// compute sum of apodization weights
				for( int index_element = index_aperture_lb; index_element <= index_aperture_ub; index_element++ )
				{

					// current lateral distance to focus
					distance_x_act = pos_rx_ctr_x[ index_element ] - pos_focus_x;

					// check location of array element
					if( distance_x_act < 0 )
					{
						// array element belongs to left aperture [ pos_rx_ctr_x[ index_element ] < pos_focus_x ]
						apodization_act = apply( window_rx, distance_x_act, width_aperture_left_over_two );
					}
					else
					{
						// array element belongs to right aperture [ pos_rx_ctr_x[ index_element ] >= pos_focus_x ]
						apodization_act = apply( window_rx, distance_x_act, width_aperture_right_over_two );
					}
					apodization_sum = apodization_sum + apodization_act;

					// store apodization weight for current rx element
					if( index_element == index_rx ) apodization = apodization_act;

				} // for( int index_element = index_aperture_lb; index_element <= index_aperture_ub; index_element++ )

				// normalize apodization weight
				if( apodization_sum > FLT_EPSILON )
				{
					apodization = apodization / apodization_sum;
				}
				else
				{
					apodization = 0.0f;
				}

			} // if( l_x < N_pos_lat_x && l_z < N_pos_lat_z )

			// argument of complex exponential
			argument_act = argument * ( index_f_lb + index_f_in_band );
			__sincosf( argument_act, &value_sine, &value_cosine );

			// assign pixel values
			pixel_value.x = pixel_value.x + apodization * ( data_RF_dft_times_two_shared[ index_f_in_block ].x * value_cosine - data_RF_dft_times_two_shared[ index_f_in_block ].y * value_sine );
			pixel_value.y = pixel_value.y + apodization * ( data_RF_dft_times_two_shared[ index_f_in_block ].x * value_sine + data_RF_dft_times_two_shared[ index_f_in_block ].y * value_cosine );

		} // for( index_f_in_block = 0; index_f_in_block < N_THREADS_PER_BLOCK; index_f_in_block++ )

		// synchronize threads to ensure complete calculations before loading new frequency data
		__syncthreads();

	} // for( index_block_f = 0; index_block_f < N_blocks_Omega_bp; index_block_f++ )

	if( l_x < N_pos_lat_x && l_z < N_pos_lat_z )
	{
		image[ l_x * N_pos_lat_z + l_z ].x += pixel_value.x;
		image[ l_x * N_pos_lat_z + l_z ].y += pixel_value.y;
	}
}