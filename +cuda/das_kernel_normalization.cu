//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// DAS kernel, normalization of monofrequent image
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// author: Martin F. Schiffner
// date: 2023-12-22
// modified: 2023-12-22

__global__ void das_kernel_normalization( t_float_complex_gpu* const image,
										  const t_float_complex_gpu* const image_f, const t_float_gpu* const weights, const t_float_gpu coefficient,
										  const int N_positions_x, const int N_positions_z )
{

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 0.) local variables
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// a) thread indices
	const int l_x = blockIdx.x * blockDim.x + threadIdx.x;
	const int l_z = blockIdx.y * blockDim.y + threadIdx.y;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// 1.) compute propagation distances
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// ensure validity of positions
	if( l_x < N_positions_x && l_z < N_positions_z )
	{
		// a) current weight
		float weight_act = weights[ l_x * N_positions_z + l_z ];

		// b) argument for complex exponential in inverse DFT
		image[ l_x * N_positions_z + l_z ].x = coefficient * image[ l_x * N_positions_z + l_z ].x + image_f[ l_x * N_positions_z + l_z ].x / weight_act;
		image[ l_x * N_positions_z + l_z ].y = coefficient * image[ l_x * N_positions_z + l_z ].y + image_f[ l_x * N_positions_z + l_z ].y / weight_act;
	}

}
