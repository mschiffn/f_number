%
% superclass for all grating lobe-derived F-numbers that enforce
% a lower bound on the angular distance
%
% author: Martin F. Schiffner
% date: 2021-08-07
% modified: 2021-09-06
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
            indicator_relevant = ( element_pitch_over_lambda > angle_lb.lower_bound ) & ( element_pitch_over_lambda < angle_lb.upper_bound );
            indicator_impossible = element_pitch_over_lambda >= angle_lb.upper_bound;

            % initialize F-number w/ zeros
            values = zeros( size( element_pitch_over_lambda ) );

            % F-number ensures lower bound on the first-order grating lobe angle
            values( indicator_relevant ) = sqrt( 1 ./ ( 1 ./ element_pitch_over_lambda( indicator_relevant ) - angle_lb.thresh ).^2 - 1 ) / 2;

            % enforce upper bound on the F-number
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
