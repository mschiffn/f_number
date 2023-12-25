//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// DAS kernel, fixed F-number
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// author: Martin F. Schiffner
// date: 2022-02-20
// modified: 2023-12-23

__global__ void das_F_number_constant( t_float_complex_gpu* const image, t_float_gpu* const weights,
									   const t_float_gpu* const distances_tx, const t_float_gpu argument_coefficient, const t_float_gpu argument_add,
									   const t_float_complex_gpu* const data_RF_dft,
									   const t_float_gpu* const positions_x, const t_float_gpu* const positions_z, const int N_positions_x, const int N_positions_z,
									   const t_float_gpu pos_rx_ctr_x, const t_float_gpu element_pitch, const t_float_gpu* const pos_rx_lb_x, const t_float_gpu* const pos_rx_ub_x, const int N_elements, const t_float_gpu M_elements,
									   const t_float_gpu f_number, t_window_ptr window_rx,
									   const int N_blocks_Omega_bp, const int index_f_lb, const int index_f_ub )
{

	// each thread block computes one sub-matrix of final DAS image

	// initialize local variables
	int l_x = blockIdx.x * blockDim.x + threadIdx.x;  //current x-index
	int l_z = blockIdx.y * blockDim.y + threadIdx.y;  //current z-index
	int index_thread = threadIdx.y * blockDim.x + threadIdx.x;

	cufftComplex pixel_value;
	pixel_value.x = 0;
	pixel_value.y = 0;

	// compute time-shift required for current x and z coordinates
	// compute argument of twiddle factor
	float argument = 0.0f;
	float apodization = 0.0f;

	float argument_act = 0.0f;
	float value_cosine = 0.0f; float value_sine = 0.0f;
	// float value_cosine_inc = 0.0f; float value_sine_inc = 0.0f;

	// shared memory
	__shared__ cufftComplex data_RF_dft_shared[ N_THREADS_PER_BLOCK ];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) compute argument of complex exponential function and apodization weight
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if( l_x < N_positions_x && l_z < N_positions_z )
	{

		//--------------------------------------------------------------------------------------------------------------------------------------
		// 1.) compute argument
		//--------------------------------------------------------------------------------------------------------------------------------------
		// a) compute round-trip distance
		float pos_focus_x = positions_x[ l_x ];
		float distance_x = pos_rx_ctr_x - pos_focus_x;
		float distance_z = positions_z[ l_z ];
		float distance_tx_rx = distances_tx[ l_x * N_positions_z + l_z ] + __fsqrt_rn( distance_x * distance_x + distance_z * distance_z );

		// b) argument for complex exponential in inverse DFT
		argument = argument_coefficient * distance_tx_rx + argument_add;

		//--------------------------------------------------------------------------------------------------------------------------------------
		// 2.) compute apodization weight
		//--------------------------------------------------------------------------------------------------------------------------------------
		// check position of selected array element
		if( distance_x < 0 )
		{
			//----------------------------------------------------------------------------------------------------------------------------------
			// selected array element belongs to left aperture [ pos_rx_ctr_x < pos_focus_x ]
			//----------------------------------------------------------------------------------------------------------------------------------
			// a) default lower bound (F-number of zero)
			int index_aperture = 0;

			// b) lower bound for positive F-number
			if( f_number > FLT_EPSILON )
			{
				float width_aperture_over_two_desired = distance_z / ( 2 * f_number );
				index_aperture = ( int ) fmaxf( ceilf( M_elements + ( pos_focus_x - width_aperture_over_two_desired ) / element_pitch ), 0 );
			}

			// c) apodization weight
			if( index_aperture < N_elements ) apodization = apply( window_rx, distance_x, pos_focus_x - pos_rx_lb_x[ index_aperture ] );
		}
		else
		{
			//----------------------------------------------------------------------------------------------------------------------------------
			// selected array element belongs to right aperture [ pos_rx_ctr_x >= pos_focus_x ]
			//----------------------------------------------------------------------------------------------------------------------------------
			// a) default upper bound (F-number of zero)
			int index_aperture = N_elements - 1;

			// b) upper bound for positive F-number
			if( f_number > FLT_EPSILON )
			{
				float width_aperture_over_two_desired = distance_z / ( 2 * f_number );
				index_aperture = ( int ) fminf( floorf( M_elements + ( pos_focus_x + width_aperture_over_two_desired ) / element_pitch ), N_elements - 1 );
			}
			
			// c) apodization weight
			if( index_aperture >= 0 ) apodization = apply( window_rx, distance_x, pos_rx_ub_x[ index_aperture ] - pos_focus_x );
		}

	} // if( l_x < N_positions_x && l_z < N_positions_z )

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 2.) loop over all frequency blocks that are required to compute sub-matrix of final image
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// __sincosf( argument, &value_sine_inc, &value_cosine_inc );
	// __sincosf( index_f_lb * argument, &value_sine, &value_cosine );

	for( int index_block_f = 0; index_block_f < N_blocks_Omega_bp; index_block_f++ )
	{

		// each thread in current block transfers one complex DFT value from device memory to shared memory

		// compute current frequency index, take into account thread ID
		int index_f_lb_block = index_f_lb + index_block_f * N_THREADS_PER_BLOCK;
		int index_f = index_f_lb_block + index_thread;

		// check validity of frequency
		if( index_f <= index_f_ub )
		{
			// transfer complex DFT value from device memory to shared memory
			// each thread loads the value for one frequency index
			data_RF_dft_shared[ index_thread ].x = data_RF_dft[ index_f ].x;
			data_RF_dft_shared[ index_thread ].y = data_RF_dft[ index_f ].y;
		}
		else
		{
			data_RF_dft_shared[ index_thread ].x = 0.0f;
			data_RF_dft_shared[ index_thread ].y = 0.0f;
		}

		// synchronize threads to ensure a complete data set before starting calculations
		__syncthreads();

		// special case: zero frequency
		if( ( index_f_lb == 0 ) && ( index_block_f == 0 ) )
		{
			pixel_value.x = - 0.5f * data_RF_dft_shared[ 0 ].x;
			pixel_value.y = - data_RF_dft_shared[ 0 ].y;
		}

		// iterate frequencies in current block
		for( int index_f_in_block = 0; index_f_in_block < N_THREADS_PER_BLOCK; index_f_in_block++ )
		{

			// argument of complex exponential
			argument_act = argument * ( index_f_lb_block + index_f_in_block );
			__sincosf( argument_act, &value_sine, &value_cosine );

			// assign pixel values
			pixel_value.x = pixel_value.x + data_RF_dft_shared[ index_f_in_block ].x * value_cosine - data_RF_dft_shared[ index_f_in_block ].y * value_sine;
			pixel_value.y = pixel_value.y + data_RF_dft_shared[ index_f_in_block ].x * value_sine + data_RF_dft_shared[ index_f_in_block ].y * value_cosine;

			// // update exponential
			// value_cosine = value_cosine * value_cosine_inc - value_sine * value_sine_inc;
			// value_sine = value_cosine * value_sine_inc + value_sine * value_cosine_inc;

		} // for( int index_f_in_block = 0; index_f_in_block < N_THREADS_PER_BLOCK; index_f_in_block++ )

		// synchronize threads to ensure complete calculations before loading new frequency data
		__syncthreads();

	} // for( index_block_f = 0; index_block_f < N_blocks_Omega_bp; index_block_f++ )

	if( l_x < N_positions_x && l_z < N_positions_z && apodization > FLT_EPSILON )
	{
		// update voxel value
		image[ l_x * N_positions_z + l_z ].x += apodization * pixel_value.x;
		image[ l_x * N_positions_z + l_z ].y += apodization * pixel_value.y;

		// update weights
		weights[ l_x * N_positions_z + l_z ] += apodization;
	}

}