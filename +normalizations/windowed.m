%
% superclass for all active window-based normalizations
%
% author: Martin F. Schiffner
% date: 2022-02-21
% modified: 2022-02-21
%
classdef windowed < normalizations.on

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)

        % independent properties
        window ( 1, 1 ) windows.window { mustBeNonempty } = windows.tukey( 0.2 )

    end % properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = windowed( windows )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % ensure class windows.window
            if ~isa( windows, 'windows.window' )
                errorStruct.message = 'windows must be windows.window!';
                errorStruct.identifier = 'windowed:NoWindow';
                error( errorStruct );
            end

            %--------------------------------------------------------------
            % 2.) create window-based normalizations
            %--------------------------------------------------------------
            % constructor of superclass
            objects@normalizations.on( size( windows ) );

            % iterate window-based normalizations
            for index_object = 1:numel( windows )

                % independent properties
                objects( index_object ).window = windows( index_object );

            end % for index_object = 1:numel( windows )

		end % function objects = windowed( windows )

    end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % apply normalization (scalar)
        %------------------------------------------------------------------
        function samples_norm = apply_scalar( windowed, samples )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class normalizations.normalization for wiener (scalar)
            % calling method ensures for samples

            %--------------------------------------------------------------
            % 2.) compute samples (scalar)
            %--------------------------------------------------------------
            % method apply of superclass
            samples_norm = apply_scalar@normalizations.on( windowed, samples );

            % apply window
            lengths_over_two = sum( samples_norm > 0, 1 ) / 2;
            temp = compute_samples( windowed.window, ( ( 0:( size( samples, 1 ) - 1 ) ).' - lengths_over_two ) ./ lengths_over_two );
            samples_norm = samples_norm .* temp;

        end % function samples_norm = apply_scalar( windowed, samples )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( windowed )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class normalizations.normalization for window (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "on_%s", windowed.window );

        end % function str_out = string_scalar( windowed )

    end % methods (Access = protected, Hidden)

end % classdef windowed < normalizations.on
