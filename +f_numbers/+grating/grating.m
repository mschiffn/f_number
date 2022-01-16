%
% abstract superclass for all grating lobe-derived F-numbers
%
% author: Martin F. Schiffner
% date: 2021-08-07
% modified: 2021-09-02
%
classdef (Abstract) grating < f_numbers.f_number

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)

        % independent properties
        angle_lb_deg ( 1, 1 ) double { mustBeNonempty, mustBeNonnegative, mustBeLessThanOrEqual( angle_lb_deg, 90 ) } = 45 % lower bound on the angular distances of the first-order grating lobes (degree)
        F_number_ub ( 1, 1 ) double { mustBeNonempty, mustBeNonnegative } = 3 % upper bound on the F-number
        anti_aliasing_tf ( 1, 1 ) logical { mustBeNonempty, mustBeNonnegative } = true % include anti-aliasing filter

        % dependent properties
        thresh ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = sin( deg2rad( 45 ) )	% threshold resulting from lower bound on the first-order grating lobe angle
        lower_bound ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 1 / ( 1 + sin( deg2rad( 45 ) ) )
        upper_bound ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 1 / ( 1 / sqrt( 1 + ( 2 * 3 )^2 ) + sin( deg2rad( 45 ) ) )
        critical ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 1 / sin( deg2rad( 45 ) )

    end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = grating( angles_lb_deg, F_numbers_ub, anti_aliasing_tfs )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at least two and at most three arguments
            narginchk( 2, 3 );

            % property validation function ensures valid angles_lb_deg
            % property validation function ensures valid F_numbers_ub

            % ensure existence of nonempty anti_aliasing_tfs
            if nargin < 3 || isempty( anti_aliasing_tfs )
                anti_aliasing_tfs = true( size( angles_lb_deg ) );
            end

            % property validation function ensures valid anti_aliasing_tfs

            % ensure equal number of dimensions and sizes
            [ angles_lb_deg, F_numbers_ub, anti_aliasing_tfs ] = auxiliary.ensureEqualSize( angles_lb_deg, F_numbers_ub, anti_aliasing_tfs );

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
                objects( index_object ).anti_aliasing_tf = anti_aliasing_tfs( index_object );

                % set dependent properties
                objects( index_object ).thresh = sin( deg2rad( objects( index_object ).angle_lb_deg ) );
%                 if objects( index_object ).anti_aliasing
%                     objects( index_object ).lower_bound = 1 / ( 2 * objects( index_object ).thresh );
%                 else
                    objects( index_object ).lower_bound = 1 / ( 1 + objects( index_object ).thresh );
%                 end
                objects( index_object ).upper_bound = 1 / ( 1 / sqrt( 1 + ( 2 * objects( index_object ).F_number_ub )^2 ) + objects( index_object ).thresh );
                objects( index_object ).critical = 1 / objects( index_object ).thresh;

            end % for index_object = 1:numel( objects )

        end % function objects = grating( angles_lb_deg, F_numbers_ub, anti_aliasing_tfs )

	end % methods

end % classdef (Abstract) grating < f_numbers.f_number
