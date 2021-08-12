%
% superclass for all boxcar (rectangular / Dirichlet) window functions
%
% author: Martin F. Schiffner
% date: 2021-08-10
% modified: 2021-08-11
%
classdef boxcar < windows.tukey

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
            % superclass ensures valid varargin

            %--------------------------------------------------------------
            % 2.) create boxcar windows
            %--------------------------------------------------------------
            % constructor of superclass
            objects@windows.tukey( zeros( size ) );

        end % function objects = boxcar( varargin )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute samples (scalar)
        %------------------------------------------------------------------
        function samples = compute_samples_scalar( ~, positions, widths_over_2 )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute samples (scalar)
            %--------------------------------------------------------------
            samples = abs( positions ) <= widths_over_2;

        end % function samples = compute_samples_scalar( ~, positions, widths_over_2 )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( ~ )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = "boxcar";

        end % function str_out = string_scalar( ~ )

	end % methods (Access = protected, Hidden)

end % classdef boxcar < windows.window
