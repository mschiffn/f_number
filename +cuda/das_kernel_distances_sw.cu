//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// DAS kernel, propagation distances of spherical wave w/ virtual source
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// author: Martin F. Schiffner
// date: 2025-08-23
// modified: 2025-08-23

__global__ void das_kernel_distances_sw( t_float_gpu* const distances,
                                         const t_float_gpu* const positions_x, const t_float_gpu* const positions_z, const int N_positions_x, const int N_positions_z,
                                         const t_float_gpu pos_tx_x, const t_float_gpu pos_tx_z )
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
		float distance_x = positions_x[ l_x ] - pos_tx_x;
		float distance_z = positions_z[ l_z ] - pos_tx_z;
		distances[ l_x * N_positions_z + l_z ] = __fsqrt_rn( distance_x * distance_x + distance_z * distance_z );
	}

}