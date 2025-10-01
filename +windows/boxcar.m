%
% superclass for all boxcar (rectangular / Dirichlet) window functions
%
% author: Martin F. Schiffner
% date: 2021-08-10
% modified: 2023-12-20
%
classdef boxcar < windows.window

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = boxcar( size )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at most one input
            narginchk( 0, 1 );

            % ensure existence of nonempty size
            if nargin < 1 || isempty( size )
                size = [ 1, 1 ];
            end

            % superclass ensures valid size

            %--------------------------------------------------------------
            % 2.) create boxcar windows
            %--------------------------------------------------------------
            % constructor of superclass
            objects@windows.window( size );

        end % function objects = boxcar( size )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute samples (scalar)
        %------------------------------------------------------------------
        function samples = compute_samples_scalar( ~, positions_over_halfwidth_abs )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window for boxcar (scalar)
            % calling method ensures positions_over_halfwidth_abs < 1

            %--------------------------------------------------------------
            % 2.) compute samples (scalar)
            %--------------------------------------------------------------
            samples = double( positions_over_halfwidth_abs < 1 );

        end % function samples = compute_samples_scalar( ~, positions_over_halfwidth_abs )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( ~ )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window for boxcar (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = "boxcar";

        end % function str_out = string_scalar( ~ )

	end % methods (Access = protected, Hidden)

end % classdef boxcar < windows.window
