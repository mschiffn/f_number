%
% superclass for all grating lobe-derived F-numbers with
% lower bound on the angle
%
% author: Martin F. Schiffner
% date: 2021-08-07
% modified: 2021-08-07
%
classdef angle_lb < f_numbers.grating.grating

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = angle_lb( varargin )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % superclass ensures valid varargin

            %--------------------------------------------------------------
            % 2.) create grating lobe-derived F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.grating.grating( varargin{ : } );

        end % function objects = angle_lb( varargin )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( angle_lb, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            % detect valid frequencies
            indicator_impossible = 1 ./ element_pitch_over_lambda <= angle_lb.thresh;

            % lower bound on the F-number (no overlap)
            F_number_lb = sqrt( max( element_pitch_over_lambda.^2, 0.25 ) - 0.25 );

            % F-number ensures lower bound on the first-order grating lobe angle
            F_number_grating = sqrt( max( 1 ./ ( 4 * ( 1 ./ element_pitch_over_lambda - angle_lb.thresh ).^2 ), 0.25 ) - 0.25 );

            % enforce lower bound on the F-number
            values = max( F_number_lb, F_number_grating );

            % enforce upper bound on the F-number
            values = min( values, angle_lb.F_number_ub );
            values( indicator_impossible ) = angle_lb.F_number_ub;

        end % function values = compute_values_scalar( angle_lb, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( angle_lb )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "angle_lb_deg_%.2f_F_ub_%.2f", angle_lb.angle_lb_deg, angle_lb.F_number_ub );

        end % function strs_out = string( angle_lb )

	end % methods (Access = protected, Hidden)

end % classdef angle_lb < f_numbers.grating.grating
