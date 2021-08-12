%
% superclass for all Tukey (cosine-tapered) windows
%
% author: Martin F. Schiffner
% date: 2021-08-10
% modified: 2021-08-11
%
classdef tukey < windows.window

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess = private)

        % independent properties
        roll_off_factor ( 1, 1 ) double { mustBeNonnegative, mustBeLessThanOrEqual( roll_off_factor, 1 ), mustBeNonempty } = 1 % roll-off factor

    end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = tukey( roll_off_factors )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % property validation functions ensure valid roll_off_factors

            %--------------------------------------------------------------
            % 2.) create Tukey windows
            %--------------------------------------------------------------
            % constructor of superclass
            objects@windows.window( size( roll_off_factors ) );

            % iterate Tukey windows
            for index_object = 1:numel( roll_off_factors )

                % set independent properties
                objects( index_object ).roll_off_factor = roll_off_factors( index_object );

            end % for index_object = 1:numel( roll_off_factors )

        end % function objects = tukey( roll_off_factors )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute samples (scalar)
        %------------------------------------------------------------------
        function samples = compute_samples_scalar( tukey, positions, widths_over_2 )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute samples (scalar)
            %--------------------------------------------------------------
            % compute lower and upper bounds
            length = tukey.roll_off_factor * widths_over_2;
            positions_abs_thresh = widths_over_2 - length;

            % absolute values of the positions
            positions_abs = abs( positions );

            % position indicators
            samples = double( positions_abs <= widths_over_2 );
            indicator_taper = ( positions_abs > positions_abs_thresh ) & samples;
            positions_abs_diff_over_length = ( positions_abs - positions_abs_thresh ) ./ length;
            samples( indicator_taper ) = ( 1 + cos( pi * positions_abs_diff_over_length( indicator_taper ) ) ) / 2;

        end % function samples = compute_values_scalar( tukey, positions, widths_over_2 )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( tukey )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "tukey_%.2f", tukey.roll_off_factor );

        end % function str_out = string_scalar( tukey )

	end % methods (Access = protected, Hidden)

end % classdef tukey < windows.window
