//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// DAS kernel, complex-valued apodization weights
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// author: Martin F. Schiffner
// date: 2023-12-14
// modified: 2025-08-26

__global__ void saft_kernel_phase_shifts( cufftComplex* const phase, t_float_gpu* const weights_tx, t_float_gpu* const weights_rx,
										  const t_float_gpu* const distances_tx, const t_float_gpu pos_tx_x, const int index_f,
										  const t_float_gpu* const positions_x, const t_float_gpu* const positions_z, const int N_positions_x, const int N_positions_z,
										  const t_float_gpu* const pos_rx_ctr_x, const t_float_gpu element_pitch, const t_float_gpu* const pos_rx_lb_x, const t_float_gpu* const pos_rx_ub_x, const int N_elements, const t_float_gpu M_elements, const int N_blocks_rx,
										  const t_float_gpu f_number_tx, t_window_ptr window_tx, const t_float_gpu window_tx_parameter,
										  const t_float_gpu f_number_rx, t_window_ptr window_rx, const t_float_gpu window_rx_parameter,
										  const t_float_gpu argument_coefficient, const t_float_gpu argument_add )
{

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 0.) local variables
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// a) initialize local variables
	const int l_x = blockIdx.x * blockDim.x + threadIdx.x;  //current x-index
	const int l_z = blockIdx.y * blockDim.y + threadIdx.y;  //current z-index
	const int index_thread = threadIdx.y * blockDim.x + threadIdx.x;  // current thread index

	float distance_tx = 0.0f;
	float value_cosine = 0.0f; float value_sine = 0.0f;
	float apodization_tx = 0.0f;
	float apodization_rx_sum = 0.0f;

	// b) shared memory
	__shared__ t_float_gpu pos_rx_ctr_x_shared[ N_THREADS_PER_BLOCK ];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) compute complex-valued apodization weights
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------------------------------------------------------
	// a) transmit propagation distance and apodization
	//------------------------------------------------------------------------------------------------------------------------------------------
	// ensure valid voxel position
	if( l_x < N_positions_x && l_z < N_positions_z )
	{
		// propagation distance of incident wave
		distance_tx = distances_tx[ l_x * N_positions_z + l_z ];

		// tx apodization
		float pos_focus_x = positions_x[ l_x ];
		float pos_focus_z = positions_z[ l_z ];
		float distance_x = pos_tx_x - pos_focus_x;

		// check position of tx array element
		if( distance_x < 0 )
		{
			//----------------------------------------------------------------------------------------------------------------------------------
			// tx array element belongs to left aperture [ pos_tx_x < pos_focus_x ]
			//----------------------------------------------------------------------------------------------------------------------------------
			// a) default lower bound (F-number of zero)
			int index_aperture = 0;

			// b) lower bound for positive F-number
			if( f_number_tx > FLT_EPSILON )
			{
				float width_aperture_over_two_desired = pos_focus_z / ( 2 * f_number_tx );
				index_aperture = ( int ) fmaxf( ceilf( M_elements + ( pos_focus_x - width_aperture_over_two_desired ) / element_pitch ), 0 );
			}

			// c) apodization weight
			if( index_aperture < N_elements ) apodization_tx = apply( window_tx, distance_x, pos_focus_x - pos_rx_lb_x[ index_aperture ], window_tx_parameter );
		}
		else
		{
			//----------------------------------------------------------------------------------------------------------------------------------
			// tx array element belongs to right aperture [ pos_tx_x >= pos_focus_x ]
			//----------------------------------------------------------------------------------------------------------------------------------
			// a) default upper bound (F-number of zero)
			int index_aperture = N_elements - 1;

			// b) upper bound for positive F-number
			if( f_number_tx > FLT_EPSILON )
			{
				float width_aperture_over_two_desired = pos_focus_z / ( 2 * f_number_tx );
				index_aperture = ( int ) fminf( floorf( M_elements + ( pos_focus_x + width_aperture_over_two_desired ) / element_pitch ), N_elements - 1 );
			}

			// c) apodization weight
			if( index_aperture >= 0 ) apodization_tx = apply( window_tx, distance_x, pos_rx_ub_x[ index_aperture ] - pos_focus_x, window_tx_parameter );
		}

	} // if( l_x < N_positions_x && l_z < N_positions_z )

	//------------------------------------------------------------------------------------------------------------------------------------------
	// b) receive apodization and propagation distances
	//------------------------------------------------------------------------------------------------------------------------------------------
	// iterate blocks of receive elements
	for( int index_block_rx = 0; index_block_rx < N_blocks_rx; index_block_rx++ )
	{

		//--------------------------------------------------------------------------------------------------------------------------------------
		// 1.) load receive element centers into shared memory
		//--------------------------------------------------------------------------------------------------------------------------------------
		// compute current receive element index
		int index_rx_lb_block = index_block_rx * N_THREADS_PER_BLOCK;
		int index_rx = index_rx_lb_block + index_thread;

		// check validity of receive element index
		if( index_rx < N_elements )
		{
			// transfer complex DFT value from device memory to shared memory
			pos_rx_ctr_x_shared[ index_thread ] = pos_rx_ctr_x[ index_rx ];
		}
		else
		{
			pos_rx_ctr_x_shared[ index_thread ] = 0.0f;
		}

		// synchronize threads to ensure a complete data set
		__syncthreads();

		//--------------------------------------------------------------------------------------------------------------------------------------
		// 2.) compute arguments
		//--------------------------------------------------------------------------------------------------------------------------------------
		// ensure valid voxel position
		if( l_x < N_positions_x && l_z < N_positions_z )
		{
			// iterate receive element indices in current block
			for( int index_rx_in_block = 0; index_rx_in_block < N_THREADS_PER_BLOCK; index_rx_in_block++ )
			{
				// compute current receive element index
				index_rx = index_rx_lb_block + index_rx_in_block;

				// ensure valid array element
				if( index_rx < N_elements )
				{
					// compute round-trip distance
					float pos_focus_x = positions_x[ l_x ];
					float distance_x = pos_rx_ctr_x_shared[ index_rx_in_block ] - pos_focus_x;
					float distance_z = positions_z[ l_z ];
					float distance_round_trip = distance_tx + __fsqrt_rn( distance_x * distance_x + distance_z * distance_z );
					float apodization_rx = 0.0f;

					// check position of current array element
					if( distance_x < 0 )
					{
						//----------------------------------------------------------------------------------------------------------------------------------
						// current array element belongs to left aperture [ pos_rx_ctr_x < pos_focus_x ]
						//----------------------------------------------------------------------------------------------------------------------------------
						// a) default lower bound (F-number of zero)
						int index_aperture = 0;

						// b) lower bound for positive F-number
						if( f_number_rx > FLT_EPSILON )
						{
							float width_aperture_over_two_desired = distance_z / ( 2 * f_number_rx );
							index_aperture = ( int ) fmaxf( ceilf( M_elements + ( pos_focus_x - width_aperture_over_two_desired ) / element_pitch ), 0 );
						}

						// c) apodization weight
						if( index_aperture < N_elements ) apodization_rx = apply( window_rx, distance_x, pos_focus_x - pos_rx_lb_x[ index_aperture ], window_rx_parameter );
					}
					else
					{
						//----------------------------------------------------------------------------------------------------------------------------------
						// current array element belongs to right aperture [ pos_rx_ctr_x >= pos_focus_x ]
						//----------------------------------------------------------------------------------------------------------------------------------
						// a) default upper bound (F-number of zero)
						int index_aperture = N_elements - 1;

						// b) upper bound for positive F-number
						if( f_number_rx > FLT_EPSILON )
						{
							float width_aperture_over_two_desired = distance_z / ( 2 * f_number_rx );
							index_aperture = ( int ) fminf( floorf( M_elements + ( pos_focus_x + width_aperture_over_two_desired ) / element_pitch ), N_elements - 1 );
						}

						// c) apodization weight
						if( index_aperture >= 0 ) apodization_rx = apply( window_rx, distance_x, pos_rx_ub_x[ index_aperture ] - pos_focus_x, window_rx_parameter );
					}

					// c) argument for complex exponential in inverse DFT
					__sincosf( index_f * ( argument_coefficient * distance_round_trip + argument_add ), &value_sine, &value_cosine );
					phase[ ( index_rx * N_positions_x + l_x ) * N_positions_z + l_z ].x = apodization_tx * apodization_rx * value_cosine;
					phase[ ( index_rx * N_positions_x + l_x ) * N_positions_z + l_z ].y = apodization_tx * apodization_rx * value_sine;

					apodization_rx_sum += apodization_rx;
				}
			}
		}

		// synchronize threads to ensure complete calculations before loading new frequency data
		__syncthreads();

	} // for( int index_block_rx = 0; index_block_rx < N_blocks_rx; index_block_rx++ )

	// ensure valid voxel position
	if( l_x < N_positions_x && l_z < N_positions_z )
	{
		weights_tx[ l_x * N_positions_z + l_z ] += apodization_tx;
		weights_rx[ l_x * N_positions_z + l_z ] = apodization_rx_sum;
	}
}