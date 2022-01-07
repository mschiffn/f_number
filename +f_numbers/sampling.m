%
% superclass for all sampling theorem-derived F-numbers with
% a fixed angular distance between
% the main lobe and
% the first-order grating lobes
%
% author: Martin F. Schiffner
% date: 2021-09-06
% modified: 2022-01-07
%
classdef sampling < f_numbers.f_number

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess = private)

        % independent properties
        distance_deg ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 5	% fixed angular distance (degree)
        F_number_ub ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 3	% upper bound on the F-number (1)

        % dependent properties
        thresh ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = sin( 5 * pi / 180 )	% additive constant resulting from the fixed angular distance (1)

    end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = sampling( distances_deg, F_numbers_ub )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure two arguments
            narginchk( 2, 2 );

            % property validation function ensures valid distances_deg
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
                objects( index_object ).thresh = sin( objects( index_object ).distance_deg * pi / 180 );

            end % for index_object = 1:numel( objects )

        end % function objects = sampling( distances_deg, F_numbers_ub )

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
            % calling method ensures class f_numbers.f_number for sampling (scalar)
            % calling method ensures nonempty positive element_pitch_over_lambda for element_pitch_over_lambda (scalar)

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            % lower bound on the F-number (no overlap)
%             F_number_lb = sqrt( max( element_pitch_over_lambda.^2, 0.25 ) - 0.25 );

            % detect valid frequencies
            indicator_possible = sampling.thresh * element_pitch_over_lambda < 1;

            % enforce upper bound on the F-number
            values = repmat( sampling.F_number_ub, size( element_pitch_over_lambda ) );

            % prevent overlap of first-order grating lobes with the main lobe
            values( indicator_possible ) = sqrt( max( 1 ./ ( 1 ./ element_pitch_over_lambda - sampling.thresh ).^2, 0.25 ) - 0.25 );

            % enforce upper bound on the F-number
            values( indicator_possible ) = min( values( indicator_possible ), sampling.F_number_ub );

        end % function values = compute_values_scalar( sampling, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( sampling )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for sampling (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "sampling_distance_deg_%.2f_F_ub_%.2f", sampling.distance_deg, sampling.F_number_ub );

        end % function str_out = string_scalar( sampling )

	end % methods (Access = protected, Hidden)

end % classdef sampling < f_numbers.f_number
