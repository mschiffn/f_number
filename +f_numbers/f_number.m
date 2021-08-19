%
% abstract superclass for all F-numbers
%
% author: Martin F. Schiffner
% date: 2021-08-03
% modified: 2021-08-19
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

            % multiple f_numbers, single element_pitch_over_lambda
            if ~isscalar( f_numbers ) && isscalar( element_pitch_over_lambda )
                element_pitch_over_lambda = repmat( element_pitch_over_lambda, size( f_numbers ) );
            end

            % single f_numbers, multiple element_pitch_over_lambda
            if isscalar( f_numbers ) && ~isscalar( element_pitch_over_lambda )
                f_numbers = repmat( f_numbers, size( element_pitch_over_lambda ) );
            end

            % ensure equal sizes
            if ~isequal( size( f_numbers ), size( element_pitch_over_lambda ) )
                errorStruct.message = 'f_numbers and element_pitch_over_lambda must have equal sizes!';
                errorStruct.identifier = 'compute_values:SizeMismatch';
                error( errorStruct );
            end

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
                % ensure positive element_pitch_over_lambda
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
        % string array (overload string method)
        %------------------------------------------------------------------
        function strs_out = string( f_numbers )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
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
