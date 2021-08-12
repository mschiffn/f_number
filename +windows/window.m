%
% abstract superclass for all window functions
%
% author: Martin F. Schiffner
% date: 2021-08-10
% modified: 2021-08-10
%
classdef (Abstract) window

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = window( size )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at most one input
            narginchk( 0, 1 );

            % ensure existence of nonempty size
            if nargin < 1 || isempty( size )
                size = [ 1, 1 ];
            end

            % ensure row vector for size
            if ~isrow( size )
                errorStruct.message = 'size must be a row vector!';
                errorStruct.identifier = 'window:NoRowVector';
                error( errorStruct );
            end

            % ensure nonempty positive integers
            mustBePositive( size );
            mustBeInteger( size );
            mustBeNonempty( size );

            %--------------------------------------------------------------
            % 2.) create windows
            %--------------------------------------------------------------
            % repeat default window
            objects = repmat( objects, size );

        end % function objects = window( size )

        %------------------------------------------------------------------
        % compute samples
        %------------------------------------------------------------------
        function samples = compute_samples( windows, positions, widths_over_2 )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure class windows.window
            if ~isa( windows, 'windows.window' )
                errorStruct.message = 'windows must be windows.window!';
                errorStruct.identifier = 'compute_samples:NoWindows';
                error( errorStruct );
            end

            % ensure cell array for positions
            if ~iscell( positions )
                positions = { positions };
            end

            % ensure cell array for widths_over_2
            if ~iscell( widths_over_2 )
                widths_over_2 = { widths_over_2 };
            end

            % TODO: ensure equal sizes
%             [ positions, widths_over_2 ] = auxiliary.mustBeEqualSizes( positions, widths_over_2 );

            %--------------------------------------------------------------
            % 2.) compute samples
            %--------------------------------------------------------------
            % specify cell array for samples
            samples = cell( size( windows ) );

            % iterate windows
            for index_object = 1:numel( windows )

                % compute samples (scalar)
                samples{ index_object } = compute_samples_scalar( windows( index_object ), positions{ index_object }, widths_over_2{ index_object } );

            end % for index_object = 1:numel( f_numbers )

            % avoid cell array for scalar windows
            if isscalar( windows )
                samples = samples{ 1 };
            end

        end % function samples = compute_samples( windows, positions, widths_over_2 )

        %------------------------------------------------------------------
        % string array (overload string method)
        %------------------------------------------------------------------
        function strs_out = string( windows )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure class windows.window
            if ~isa( windows, 'windows.window' )
                errorStruct.message = 'windows must be windows.window!';
                errorStruct.identifier = 'string:NoWindows';
                error( errorStruct );
            end

            %--------------------------------------------------------------
            % 2.) string array
            %--------------------------------------------------------------
            % repeat empty string
            strs_out = repmat( "", size( windows ) );

            % iterate windows
            for index_object = 1:numel( windows )

                strs_out( index_object ) = string_scalar( windows( index_object ) );

            end % for index_object = 1:numel( windows )

        end % function strs_out = string( windows )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (Abstract, protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Abstract, Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute samples (scalar)
        %------------------------------------------------------------------
        samples = compute_samples_scalar( window, positions, widths_over_2 );

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        str_out = string_scalar( window );

	end % methods (Abstract, Access = protected, Hidden)

end % classdef (Abstract) window
