%
% superclass for all sampling theorem-derived F-numbers with
% an upper bound on the angle
%
% author: Martin F. Schiffner
% date: 2021-09-06
% modified: 2021-09-06
%
classdef sampling < f_numbers.f_number

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess = private)

        % independent properties
        angle_ub_deg ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 10     % upper bound on the angle (degree)
        F_number_ub ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 3       % upper bound on the F-number

        % dependent properties
        thresh ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 2 * tan( 10 * pi / 180 )	% threshold resulting from lower bound on the first-order grating lobe angle

	end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = sampling( angles_lb_deg, F_numbers_ub )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % superclass ensures valid varargin

            %--------------------------------------------------------------
            % 2.) create sampling theorem-derived F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.f_number( size( angles_lb_deg ) );

            % iterate sampling theorem-derived F-numbers
            for index_object = 1:numel( objects )

                % set independent properties
                objects( index_object ).angle_ub_deg = angles_lb_deg( index_object );
                objects( index_object ).F_number_ub = F_numbers_ub( index_object );

                % set dependent properties
                objects( index_object ).thresh = 2 * tan( objects( index_object ).angle_ub_deg * pi / 180 );

            end % for index_object = 1:numel( objects )

        end % function objects = sampling( angles_lb_deg, F_numbers_ub )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( sampling, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            % lower bound on the F-number (no overlap)
            F_number_lb = sqrt( max( element_pitch_over_lambda.^2, 0.25 ) - 0.25 );

            % detect valid frequencies
            indicator_impossible = sampling.angle_ub_deg * pi / 180 >= atan( 1 ./ ( 2 * F_number_lb ) );

            % F-number ensures lower bound on the first-order grating lobe angle
            values = 1 ./ ( 1./ F_number_lb - sampling.thresh );

            % enforce upper bound on the F-number
            values = min( values, sampling.F_number_ub );
            values( indicator_impossible ) = sampling.F_number_ub;

        end % function values = compute_values_scalar( sampling, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( sampling )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "sampling_deg_%.2f_F_ub_%.2f", sampling.angle_ub_deg, sampling.F_number_ub );

        end % function strs_out = string( sampling )

	end % methods (Access = protected, Hidden)

end % classdef sampling < f_numbers.f_number
