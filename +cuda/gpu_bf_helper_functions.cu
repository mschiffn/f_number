//-------------------------------------------------------------------------
// helper functions
//-------------------------------------------------------------------------
// #include "gpu_bf_error_handling.cuh"

//-------------------------------------------------------------------------
// call "isa" in MATLAB workspace to make sure subclasses are identified
//-------------------------------------------------------------------------
// workaround to check MATLAB classes
bool isa( const mxArray* object, const char* type )
{
	//Create LHS/RHS arrays for calling MATLAB
	mxArray* lhs[ 1 ];
	mxArray* rhs[ 2 ];
	//Return value
	bool retVal;
	//Populate Inputs to MATLAB isa
	rhs[ 0 ] = mxDuplicateArray( object );	// make deep copy of object to avoid problems w/ const modifier
	rhs[ 1 ] = mxCreateString( type );
	//Call the MATLAB isa function
	mexCallMATLAB( 1, lhs, 2, rhs, "isa" );
	//Extract result
	retVal = mxIsLogicalScalarTrue( lhs[ 0 ] );
	//Cleanup
	mxDestroyArray( lhs[ 0 ] );
	mxDestroyArray( rhs[ 0 ] );
	mxDestroyArray( rhs[ 1 ] );
	//Done
	return retVal;
}

//-------------------------------------------------------------------------
// compute F-numbers
//-------------------------------------------------------------------------
double* compute_values( const mxArray* f_number, mxArray* element_pitch_norm )
{
	mxArray* lhs[ 1 ];
	mxArray* rhs[ 2 ];
	rhs[ 0 ] = mxDuplicateArray( f_number );	// make deep copy of F-number to avoid problems w/ const modifier
	rhs[ 1 ] = element_pitch_norm;

	// call method compute_values in MATLAB
	mexCallMATLAB( 1, lhs, 2, rhs, "compute_values" );

	// clean up
	mxDestroyArray( rhs[ 0 ] );
	mxDestroyArray( rhs[ 1 ] );

	// get values
	return ( double * ) mxGetData( lhs[ 0 ] );
}

//-------------------------------------------------------------------------
// map MATLAB window classes to ids
//-------------------------------------------------------------------------
int get_window_id( const mxArray* window )
{

	if( mxIsClass( window, "windows.boxcar" ) )
	{
		return 0;
	}
	else if( mxIsClass( window, "windows.hann" ) )
	{
		return 1;
	}
	else if( mxIsClass( window, "windows.tukey" ) )
	{
		// ensure existence of property fraction_cosine
		if( mxGetProperty( window, 0, "fraction_cosine" ) != NULL )
		{
			// read property fraction_cosine
			if( *( ( double * ) mxGetData( mxGetProperty( window, 0, "fraction_cosine" ) ) ) == 0.2 )
			{
				return 2;
			}
		}
	}
	else if( mxIsClass( window, "windows.triangular" ) )
	{
		return 4;
	}

	return -1;

} // int get_window_id( const mxArray* window )

//-------------------------------------------------------------------------
// convert object to string
//-------------------------------------------------------------------------
// Matlab's String class is encapsulated,
// use Matlab call to convert it to char array
//   mxArray *string_class[1], *char_array[1];
//   string_class[0] = pr;
//   mexCallMATLAB(1, char_array, 1, string_class, "char");
// Parse the char array to create an std::string
//   int buflen = mxGetN(char_array[0])*sizeof(mxChar)+1;
//   char* buf = new char[buflen];
//   mxGetString(char_array[0],buf,buflen);
//   data = std::string(buf);
//   delete buf;
int string( const mxArray* object, char* string, mwSize strlen )
{
	// if( !mxChar( object ) )
	// local variables
	mxArray* lhs_1;
	mxArray* lhs_2;
	mxArray* rhs_1;
	rhs_1 = mxDuplicateArray( object );	// make deep copy of F-number to avoid problems w/ const modifier

	// call method string in MATLAB
	mexCallMATLAB( 1, &lhs_1, 1, &rhs_1, "string" );

	// call method char in MATLAB
	mexCallMATLAB( 1, &lhs_2, 1, &lhs_1, "char" );

	// clean up
	mxDestroyArray( rhs_1 );
	mxDestroyArray( lhs_1 );

	// check for mxChar array
	// if( !mxIsChar( lhs_2 ) ) return 1;

	// mxChar array to C-style string
	return mxGetString( lhs_2, string, strlen );

}
