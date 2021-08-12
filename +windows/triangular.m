%
% superclass for all triangular (Bartlett / Fejér) windows
%
% author: Martin F. Schiffner
% date: 2021-08-11
% modified: 2021-08-11
%
classdef triangular < windows.window

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = triangular( varargin )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % property validation functions ensure valid roll_off_factors

            %--------------------------------------------------------------
            % 2.) create triangular windows
            %--------------------------------------------------------------
            % constructor of superclass
            objects@windows.window( varargin{ : } );

        end % function objects = triangular( roll_off_factors )

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
            % absolute values of the positions
            positions_abs = abs( positions );

            % compute samples
            samples = ( positions_abs <= widths_over_2 ) .* ( 1 - positions_abs ./ widths_over_2 );

        end % function samples = compute_values_scalar( ~, positions, widths_over_2 )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( triangular )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = "triangular";

        end % function str_out = string_scalar( triangular )

	end % methods (Access = protected, Hidden)

end % classdef triangular < windows.window
