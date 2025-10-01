%
% abstract superclass for all window functions
%
% author: Martin F. Schiffner
% date: 2021-08-10
% modified: 2023-12-20
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
            % calling method ensures existence of nonempty size

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
        function samples = compute_samples( windows, positions_over_halfwidth )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure two arguments
            narginchk( 2, 2 );

            % ensure class windows.window
            if ~isa( windows, 'windows.window' )
                errorStruct.message = 'windows must be windows.window!';
                errorStruct.identifier = 'compute_samples:NoWindows';
                error( errorStruct );
            end

            % ensure cell array for positions_over_halfwidth
            if ~iscell( positions_over_halfwidth )
                positions_over_halfwidth = { positions_over_halfwidth };
            end

            % ensure equal number of dimensions and sizes
            [ windows, positions_over_halfwidth ] = auxiliary.ensureEqualSize( windows, positions_over_halfwidth );

            %--------------------------------------------------------------
            % 2.) compute samples
            %--------------------------------------------------------------
            % specify cell array for samples
            samples = cell( size( windows ) );

            % iterate windows
            for index_object = 1:numel( windows )

                % absolute values of the positions
                positions_over_halfwidth_abs = abs( positions_over_halfwidth{ index_object } );

                % compute samples (scalar)
                samples{ index_object } = compute_samples_scalar( windows( index_object ), positions_over_halfwidth_abs );

            end % for index_object = 1:numel( f_numbers )

            % avoid cell array for scalar windows
            if isscalar( windows )
                samples = samples{ 1 };
            end

        end % function samples = compute_samples( windows, positions_over_halfwidth )

        %------------------------------------------------------------------
        % TODO: compute derivatives
        %------------------------------------------------------------------
        

        %------------------------------------------------------------------
        % string array (overload string method)
        %------------------------------------------------------------------
        function strs_out = string( windows )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

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
        samples = compute_samples_scalar( window, positions_over_halfwidth_abs );

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        str_out = string_scalar( window );

	end % methods (Abstract, Access = protected, Hidden)

end % classdef (Abstract) window
