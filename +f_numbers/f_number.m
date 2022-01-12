%
% abstract superclass for all F-numbers
%
% author: Martin F. Schiffner
% date: 2021-08-03
% modified: 2022-01-12
%
classdef (Abstract) f_number

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = f_number( size )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % ensure row vector for size
            if ~isrow( size )
                errorStruct.message = 'size must be a row vector!';
                errorStruct.identifier = 'f_number:NoRowVector';
                error( errorStruct );
            end

            % ensure nonempty positive integers
            mustBePositive( size );
            mustBeInteger( size );
            mustBeNonempty( size );

            %--------------------------------------------------------------
            % 2.) create F-numbers
            %--------------------------------------------------------------
            % repeat default F-number
            objects = repmat( objects, size );

        end % function objects = f_number( size )

        %------------------------------------------------------------------
        % compute values
        %------------------------------------------------------------------
        function values = compute_values( f_numbers, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure two arguments
            narginchk( 2, 2 );

            % ensure class f_numbers.f_number
            if ~isa( f_numbers, 'f_numbers.f_number' )
                errorStruct.message = 'f_numbers must be f_numbers.f_number!';
                errorStruct.identifier = 'compute_values:NoFNumbers';
                error( errorStruct );
            end

            % ensure cell array for element_pitch_over_lambda
            if ~iscell( element_pitch_over_lambda )
                element_pitch_over_lambda = { element_pitch_over_lambda };
            end

            % ensure equal number of dimensions and sizes
            [ f_numbers, element_pitch_over_lambda ] = auxiliary.ensureEqualSize( f_numbers, element_pitch_over_lambda );

            %--------------------------------------------------------------
            % 2.) compute values
            %--------------------------------------------------------------
            % specify cell array for values
            values = cell( size( f_numbers ) );

            % iterate F-numbers
            for index_object = 1:numel( f_numbers )

                %----------------------------------------------------------
                % a) check arguments
                %----------------------------------------------------------
                % ensure nonempty positive element_pitch_over_lambda{ index_object }
                mustBeNonempty( element_pitch_over_lambda{ index_object } );
                mustBePositive( element_pitch_over_lambda{ index_object } );

                %----------------------------------------------------------
                % b) compute values (scalar)
                %----------------------------------------------------------
                values{ index_object } = compute_values_scalar( f_numbers( index_object ), element_pitch_over_lambda{ index_object } );

            end % for index_object = 1:numel( f_numbers )

            % avoid cell array for single f_numbers
            if isscalar( f_numbers )
                values = values{ 1 };
            end

        end % function values = compute_values( f_numbers, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % compute bounds on the main lobe
        %------------------------------------------------------------------
        function bounds = compute_bounds_main( f_numbers, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure two arguments
            narginchk( 2, 2 );

            % method compute_values ensures class f_numbers.f_number for f_numbers
            % method compute_values ensures cell array for element_pitch_over_lambda

            % method compute_values ensures equal number of dimensions and sizes

            %--------------------------------------------------------------
            % 2.) compute bounds on the main lobe
            %--------------------------------------------------------------
            % compute values
            values = compute_values( f_numbers, element_pitch_over_lambda );

            % ensure cell array for values
            if ~iscell( values )
                values = { values };
            end

            % specify cell array for bounds
            bounds = cell( size( f_numbers ) );

            % iterate F-numbers
            for index_object = 1:numel( f_numbers )

                %----------------------------------------------------------
                % a) check arguments
                %----------------------------------------------------------
                % function compute values ensured nonempty positive element_pitch_over_lambda{ index_object }

                %----------------------------------------------------------
                % b) compute bounds on the main lobe (scalar)
                %----------------------------------------------------------
                bounds{ index_object } = 1 ./ sqrt( 1 + ( 2 * values{ index_object } ).^2 );

            end % for index_object = 1:numel( f_numbers )

            % avoid cell array for single f_numbers
            if isscalar( f_numbers )
                bounds = bounds{ 1 };
            end

        end % function bounds = compute_bounds_main( f_numbers, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % compute distances of the first-order grating lobes
        %------------------------------------------------------------------
        function distances = compute_distances_grating( f_numbers, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure two arguments
            narginchk( 2, 2 );

            % method compute_bounds_main ensures class f_numbers.f_number for f_numbers

            % ensure cell array for element_pitch_over_lambda
            if ~iscell( element_pitch_over_lambda )
                element_pitch_over_lambda = { element_pitch_over_lambda };
            end

            % ensure equal number of dimensions and sizes
            [ f_numbers, element_pitch_over_lambda ] = auxiliary.ensureEqualSize( f_numbers, element_pitch_over_lambda );

            %--------------------------------------------------------------
            % 2.) compute distances of the first-order grating lobes
            %--------------------------------------------------------------
            % compute values
            bounds = compute_bounds_main( f_numbers, element_pitch_over_lambda );

            % ensure cell array for bounds
            if ~iscell( bounds )
                bounds = { bounds };
            end

            % specify cell array for distances
            distances = cell( size( f_numbers ) );

            % iterate F-numbers
            for index_object = 1:numel( f_numbers )

                %----------------------------------------------------------
                % a) check arguments
                %----------------------------------------------------------
                % function compute values ensured nonempty positive element_pitch_over_lambda{ index_object }

                %----------------------------------------------------------
                % b) compute distances of the first-order grating lobes (scalar)
                %----------------------------------------------------------
                distances{ index_object } = 1 ./ element_pitch_over_lambda{ index_object } - bounds{ index_object };

            end % for index_object = 1:numel( f_numbers )

            % avoid cell array for single f_numbers
            if isscalar( f_numbers )
                distances = distances{ 1 };
            end

        end % function distances = compute_distances_grating( f_numbers, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (overload string method)
        %------------------------------------------------------------------
        function strs_out = string( f_numbers )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % ensure class f_numbers.f_number
            if ~isa( f_numbers, 'f_numbers.f_number' )
                errorStruct.message = 'f_numbers must be f_numbers.f_number!';
                errorStruct.identifier = 'string:NoFNumbers';
                error( errorStruct );
            end

            %--------------------------------------------------------------
            % 2.) string array
            %--------------------------------------------------------------
            % repeat empty string
            strs_out = repmat( "", size( f_numbers ) );

            % iterate F-numbers
            for index_object = 1:numel( f_numbers )

                strs_out( index_object ) = string_scalar( f_numbers( index_object ) );

            end % for index_object = 1:numel( f_numbers )

        end % function strs_out = string( f_numbers )

    end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (Abstract, protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Abstract, Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        samples = compute_values_scalar( f_number, element_pitch_over_lambda );

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        str_out = string_scalar( f_number );

	end % methods (Abstract, Access = protected, Hidden)

end % classdef (Abstract) f_number
