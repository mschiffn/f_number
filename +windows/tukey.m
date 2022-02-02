%
% superclass for all Tukey (cosine-tapered) windows
%
% author: Martin F. Schiffner
% date: 2021-08-10
% modified: 2022-02-02
%
classdef tukey < windows.window

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess = private)

        % independent properties
        roll_off_factor ( 1, 1 ) double { mustBePositive, mustBeLessThan( roll_off_factor, 1 ), mustBeNonempty } = 0.5 % roll-off factor

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
            % ensure at most one argument
            narginchk( 0, 1 );

            % ensure existence of nonempty roll_off_factors
            if nargin < 1 || isempty( roll_off_factors )
                roll_off_factors = 0.5;
            end

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
        function samples = compute_samples_scalar( tukey, positions_over_halfwidth )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for tukey (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute samples (scalar)
            %--------------------------------------------------------------
            % compute lower and upper bounds
            positions_over_halfwidth_abs_thresh = 1 - tukey.roll_off_factor;

            % absolute values of the positions
            positions_over_halfwidth_abs = abs( positions_over_halfwidth );

            % position indicators
            samples = double( positions_over_halfwidth_abs <= 1 );
            indicator_taper = ( positions_over_halfwidth_abs > positions_over_halfwidth_abs_thresh ) & samples;
            positions_over_halfwidth_abs_diff_over_length = ( positions_over_halfwidth_abs - positions_over_halfwidth_abs_thresh ) ./ tukey.roll_off_factor;
            samples( indicator_taper ) = ( 1 + cos( pi * positions_over_halfwidth_abs_diff_over_length( indicator_taper ) ) ) / 2;

        end % function samples = compute_samples_scalar( tukey, positions_over_halfwidth )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( tukey )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window for tukey (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "tukey_%.2f", tukey.roll_off_factor );

        end % function str_out = string_scalar( tukey )

	end % methods (Access = protected, Hidden)

end % classdef tukey < windows.window
