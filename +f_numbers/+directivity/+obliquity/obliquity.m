%
% abstract superclass for all obliquity factors
%
% author: Martin F. Schiffner
% date: 2022-01-08
% modified: 2022-01-08
%
classdef (Abstract) obliquity

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = obliquity( size )

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
            % 2.) create obliquity factors
            %--------------------------------------------------------------
            % repeat default obliquity factor
            objects = repmat( objects, size );

        end % function objects = obliquity( size )

        %------------------------------------------------------------------
        % compute values
        %------------------------------------------------------------------
        function values = compute_values( obliquities, thetas )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure two arguments
            narginchk( 2, 2 );

            % ensure class f_numbers.directivity.obliquity.obliquity
            if ~isa( obliquities, 'f_numbers.directivity.obliquity.obliquity' )
                errorStruct.message = 'obliquities must be f_numbers.directivity.obliquity.obliquity!';
                errorStruct.identifier = 'compute_values:NoObliquityFactors';
                error( errorStruct );
            end

            % ensure cell array for thetas
            if ~iscell( thetas )
                thetas = { thetas };
            end

            % ensure equal number of dimensions and sizes
            [ obliquities, thetas ] = auxiliary.ensureEqualSize( obliquities, thetas );

            %--------------------------------------------------------------
            % 2.) compute values
            %--------------------------------------------------------------
            % specify cell array for values
            values = cell( size( obliquities ) );

            % iterate obliquity factors
            for index_object = 1:numel( obliquities )

                %----------------------------------------------------------
                % a) check arguments
                %----------------------------------------------------------
                % ensure nonempty thetas{ index_object }
                mustBeNonempty( thetas{ index_object } );
                mustBeInRange( thetas{ index_object }, 0, pi / 2 );

                %----------------------------------------------------------
                % b) compute values (scalar)
                %----------------------------------------------------------
                values{ index_object } = compute_values_scalar( obliquities( index_object ), thetas{ index_object } );

            end % for index_object = 1:numel( obliquities )

            % avoid cell array for single obliquity factor
            if isscalar( obliquities )
                values = values{ 1 };
            end

        end % function values = compute_values( obliquities, thetas )

        %------------------------------------------------------------------
        % string array (overload string method)
        %------------------------------------------------------------------
        function strs_out = string( obliquities )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % ensure class f_numbers.directivity.obliquity.obliquity
            if ~isa( obliquities, 'f_numbers.directivity.obliquity.obliquity' )
                errorStruct.message = 'obliquities must be f_numbers.directivity.obliquity.obliquity!';
                errorStruct.identifier = 'compute_values:NoObliquityFactors';
                error( errorStruct );
            end

            %--------------------------------------------------------------
            % 2.) string array
            %--------------------------------------------------------------
            % repeat empty string
            strs_out = repmat( "", size( obliquities ) );

            % iterate obliquity factors
            for index_object = 1:numel( obliquities )

                strs_out( index_object ) = string_scalar( obliquities( index_object ) );

            end % for index_object = 1:numel( f_numbers )

        end % function strs_out = string( obliquities )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (Abstract, protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Abstract, Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        samples = compute_values_scalar( obliquity, theta );

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        str_out = string_scalar( obliquity );

	end % methods (Abstract, Access = protected, Hidden)

end % classdef (Abstract) obliquity
