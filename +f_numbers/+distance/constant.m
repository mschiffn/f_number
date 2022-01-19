%
% superclass for all F-numbers that maintain
% a fixed angular distance between
% the main lobe and
% the first-order grating lobes
%
% author: Martin F. Schiffner
% date: 2021-09-06
% modified: 2022-01-19
%
classdef constant < f_numbers.f_number

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)

        % independent properties
        distance_deg ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty, mustBeLessThanOrEqual( distance_deg, 90 ) } = 5	% fixed angular distance (degree)
        F_number_ub ( 1, 1 ) double { mustBePositive, mustBeNonempty } = 3	% upper bound on the F-number (1)

        % dependent properties
        distance_rad ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = deg2rad( 5 )                       % fixed angular distance (rad)
        sin_distance ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = sin( deg2rad( 5 ) )                % additive constant resulting from the fixed angular distance (1)
        one_plus_cos_distance ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 1 + cos( deg2rad( 5 ) )	% additive constant resulting from the fixed angular distance (1)

    end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = constant( distances_deg, F_numbers_ub )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at least one and at most two arguments
            narginchk( 1, 2 );

            % property validation function ensures valid distances_deg

            % ensure existence of nonempty F_numbers_ub
            if nargin < 2 || isempty( F_numbers_ub )
                F_numbers_ub = inf;
            end

            % property validation function ensures valid F_numbers_ub

            % ensure equal number of dimensions and sizes
            [ distances_deg, F_numbers_ub ] = auxiliary.ensureEqualSize( distances_deg, F_numbers_ub );

            %--------------------------------------------------------------
            % 2.) create sampling theorem-derived F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.f_number( size( distances_deg ) );

            % iterate sampling theorem-derived F-numbers
            for index_object = 1:numel( objects )

                % set independent properties
                objects( index_object ).distance_deg = distances_deg( index_object );
                objects( index_object ).F_number_ub = F_numbers_ub( index_object );

                % set dependent properties
                objects( index_object ).distance_rad = deg2rad( objects( index_object ).distance_deg );
                objects( index_object ).sin_distance = sin( objects( index_object ).distance_rad );
                objects( index_object ).one_plus_cos_distance = 1 + cos( objects( index_object ).distance_rad );

            end % for index_object = 1:numel( objects )

        end % function objects = constant( distances_deg, F_numbers_ub )

	end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods (protected, hidden)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( constant, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for constant (scalar)
            % calling method ensures nonempty positive element_pitch_over_lambda for element_pitch_over_lambda (scalar)

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            % detect valid frequencies
            indicator_lb = element_pitch_over_lambda >= 1 / constant.one_plus_cos_distance;
            indicator_ub = element_pitch_over_lambda < 1 / constant.sin_distance;
            indicator_transition = ~indicator_lb & element_pitch_over_lambda >= 0.5;
            indicator_possible = indicator_lb & indicator_ub;

            % initialize F-number w/ zeros
            values = zeros( size( element_pitch_over_lambda ) );

            % maintain 90° angular distance of the first-order grating lobes
            values( indicator_transition ) = sqrt( 1 ./ ( 1 ./ element_pitch_over_lambda( indicator_transition ) - 1 ).^2 - 1 ) / 2;

            % maintain fixed angular distance between the main lobe and the first-order grating lobes
            temp = element_pitch_over_lambda( indicator_possible ).^2;
            values( indicator_possible ) = ( sqrt( 2 * constant.one_plus_cos_distance * temp - 1 ) + constant.one_plus_cos_distance * constant.sin_distance * temp ) ./ ( 2 - 2 * constant.sin_distance^2 * temp );

            % enforce upper bound on the F-number
            values( indicator_possible ) = min( values( indicator_possible ), constant.F_number_ub );
            values( ~indicator_ub ) = constant.F_number_ub;

        end % function values = compute_values_scalar( constant, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( constant )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for constant (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "distance_deg_%.2f_F_ub_%.2f", constant.distance_deg, constant.F_number_ub );

        end % function str_out = string_scalar( constant )

	end % methods (Access = protected, Hidden)

end % classdef constant < f_numbers.f_number
