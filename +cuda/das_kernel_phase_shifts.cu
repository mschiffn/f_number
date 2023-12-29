//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// DAS kernel, complex-valued apodization weights
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// author: Martin F. Schiffner
// date: 2023-12-14
// modified: 2023-12-28

__global__ void das_kernel_phase_shifts( cufftComplex* const phase, t_float_gpu* const weights, const t_float_gpu* const distances_tx, const int index_f,
									  const t_float_gpu* const positions_x, const t_float_gpu* const positions_z, const int N_positions_x, const int N_positions_z,
									  const t_float_gpu* const pos_rx_ctr_x, const t_float_gpu element_pitch, const t_float_gpu* const pos_rx_lb_x, const t_float_gpu* const pos_rx_ub_x, const int N_elements, const t_float_gpu M_elements, const int N_blocks_rx,
									  const t_float_gpu f_number, t_window_ptr window_rx, const t_float_gpu window_rx_parameter,
									  const t_float_gpu argument_coefficient, const t_float_gpu argument_add )
{

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 0.) local variables
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// a) initialize local variables
	const int l_x = blockIdx.x * blockDim.x + threadIdx.x;  //current x-index
	const int l_z = blockIdx.y * blockDim.y + threadIdx.y;  //current z-index
	const int index_thread = threadIdx.y * blockDim.x + threadIdx.x;
	float distances_tx_act = 0.0f;
	float value_cosine = 0.0f; float value_sine = 0.0f;
	float apodization_sum = 0.0f;

	// b) shared memory
	__shared__ t_float_gpu pos_rx_ctr_x_shared[ N_THREADS_PER_BLOCK ];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) compute arguments
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// //--------------------------------------------------------------------------------------------------------------------------------------
	// // a) load propagation distances of incident wave into shared memory
	// //--------------------------------------------------------------------------------------------------------------------------------------
	// ensure validity of positions
	if( l_x < N_positions_x && l_z < N_positions_z )
	{
		distances_tx_act = distances_tx[ l_x * N_positions_z + l_z ];
	}

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
		// ensure validity of positions
		if( l_x < N_positions_x && l_z < N_positions_z )
		{
			// iterate receive element indices in current block
			for( int index_rx_in_block = 0; index_rx_in_block < N_THREADS_PER_BLOCK; index_rx_in_block++ )
			{

				// a) compute current receive element index
				index_rx = index_rx_lb_block + index_rx_in_block;

				if( index_rx < N_elements )
				{
					// a) compute round-trip distance
					float pos_focus_x = positions_x[ l_x ];
					float distance_x = pos_rx_ctr_x_shared[ index_rx_in_block ] - pos_focus_x;
					float distance_z = positions_z[ l_z ];
					float distance_round_trip = distances_tx_act + __fsqrt_rn( distance_x * distance_x + distance_z * distance_z );
					float apodization = 0.0f;

					// check position of selected array element
					if( distance_x < 0 )
					{
						//----------------------------------------------------------------------------------------------------------------------------------
						// current array element belongs to left aperture [ pos_rx_ctr_x < pos_focus_x ]
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
						if( index_aperture < N_elements ) apodization = apply( window_rx, distance_x, pos_focus_x - pos_rx_lb_x[ index_aperture ], window_rx_parameter );
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
							float width_aperture_over_two_desired = positions_z[ l_z ] / ( 2 * f_number );
							index_aperture = ( int ) fminf( floorf( M_elements + ( pos_focus_x + width_aperture_over_two_desired ) / element_pitch ), N_elements - 1 );
						}

						// c) apodization weight
						if( index_aperture >= 0 ) apodization = apply( window_rx, distance_x, pos_rx_ub_x[ index_aperture ] - pos_focus_x, window_rx_parameter );
					}

					// c) argument for complex exponential in inverse DFT
					__sincosf( index_f * ( argument_coefficient * distance_round_trip + argument_add ), &value_sine, &value_cosine );
					phase[ ( index_rx * N_positions_x + l_x ) * N_positions_z + l_z ].x = apodization * value_cosine;
					phase[ ( index_rx * N_positions_x + l_x ) * N_positions_z + l_z ].y = apodization * value_sine;

					apodization_sum += apodization;
				}
			}
		}

		// synchronize threads to ensure complete calculations before loading new frequency data
		__syncthreads();

	} // for( int index_block_rx = 0; index_block_rx < N_blocks_rx; index_block_rx++ )

	// ensure validity of positions
	if( l_x < N_positions_x && l_z < N_positions_z )
	{
		weights[ l_x * N_positions_z + l_z ] = apodization_sum;
	}

}