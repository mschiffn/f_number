%
% abstract superclass for all grating lobe-derived F-numbers
%
% author: Martin F. Schiffner
% date: 2021-08-07
% modified: 2021-08-07
%
classdef (Abstract) grating < f_numbers.f_number

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess = private)

        % independent properties
        angle_lb_deg ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 40             % lower bound on the first-order grating lobe angle (degree)
        F_number_ub ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 3               % upper bound on the F-number (why?)

        % dependent properties
        thresh ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = sin( 40 * pi / 180 )	% threshold resulting from lower bound on the first-order grating lobe angle

	end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = grating( angles_lb_deg, F_numbers_ub )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 2, 2 );

            % property validation functions ensure valid angles_lb

            %--------------------------------------------------------------
            % 2.) create grating lobe-derived F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.f_number( size( angles_lb_deg ) );

            % iterate grating lobe-derived F-numbers
            for index_object = 1:numel( objects )

                % set independent properties
                objects( index_object ).angle_lb_deg = angles_lb_deg( index_object );
                objects( index_object ).F_number_ub = F_numbers_ub( index_object );

                % set dependent properties
                objects( index_object ).thresh = sin( objects( index_object ).angle_lb_deg * pi / 180 );

            end % for index_object = 1:numel( objects )

        end % function objects = grating( angles_lb_deg, F_numbers_ub )

	end % methods

end % classdef (Abstract) grating < f_numbers.f_number
