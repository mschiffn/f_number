%
% superclass for all directivity-derived F-numbers
%
% author: Martin F. Schiffner
% date: 2021-08-06
% modified: 2022-01-08
%
classdef directivity < f_numbers.f_number

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)

        % independent properties
        width_over_pitch ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty, mustBeLessThanOrEqual( width_over_pitch, 1 ) } = 0.918	% element width-to-element pitch ratio (filling factor)
        obliquity_factor ( 1, 1 ) f_numbers.directivity.obliquity.obliquity { mustBeNonempty } = f_numbers.directivity.obliquity.soft	% obliquity factor for selected boundary condition
        attenuation_dB ( 1, 1 ) double { mustBePositive, mustBeNonempty } = 3                                                           % maximum permissible attenuation (dB)

        % dependent properties
        thresh ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 10^( -3 / 20 )	% threshold resulting from maximum permissible attenuation

    end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = directivity( widths_over_pitch, obliquity_factors, attenuations_dB )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at least one and at most three arguments
            narginchk( 1, 3 );

            % property validation function ensures valid widths_over_pitch

            % ensure existence of nonempty obliquity_factors
            if nargin < 2 || isempty( obliquity_factors )
                obliquity_factors = f_numbers.directivity.obliquity.soft;
            end

            % property validation function ensures valid obliquity_factors

            % ensure existence of nonempty attenuations_dB
            if nargin < 3 || isempty( attenuations_dB )
                attenuations_dB = 3;
            end

            % property validation function ensures valid attenuations_dB

            % ensure equal number of dimensions and sizes
            [ widths_over_pitch, obliquity_factors, attenuations_dB ] = auxiliary.ensureEqualSize( widths_over_pitch, obliquity_factors, attenuations_dB );

            %--------------------------------------------------------------
            % 2.) create directivity-derived F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.f_number( size( widths_over_pitch ) );

            % iterate directivity-derived F-numbers
            for index_object = 1:numel( objects )

                % set independent properties
                objects( index_object ).width_over_pitch = widths_over_pitch( index_object );
                objects( index_object ).obliquity_factor = obliquity_factors( index_object );
                objects( index_object ).attenuation_dB = attenuations_dB( index_object );

                % set dependent properties
                objects( index_object ).thresh = 10^( -objects( index_object ).attenuation_dB / 20 );

            end % for index_object = 1:numel( objects )

        end % function objects = directivity( widths_over_pitch, obliquity_factors, attenuations_dB )

	end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods (protected, hidden)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( directivity, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for directivity (scalar)
            % calling method ensures positive element_pitch_over_lambda for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            % normalized element width
            element_width_over_lambda = directivity.width_over_pitch * element_pitch_over_lambda;

            % directivity pattern
            directivity_pattern = @( theta, width_norm ) compute_values( directivity.obliquity_factor, theta ) * sinc( width_norm * sin( theta ) );

            % function to minimize
            f = @( theta, width_norm ) abs( directivity_pattern( theta, width_norm ) - directivity.thresh );

            % initialize F-numbers w/ zeros
            values = zeros( size( element_pitch_over_lambda ) );

            % iterate samples
            for index_width = 1:numel( element_pitch_over_lambda )
                alpha = fminbnd( @( theta ) f( theta, element_width_over_lambda( index_width ) ), 0, pi / 2 );
                values( index_width ) = 1 / ( 2 * tan( alpha ) );
            end

        end % function values = compute_values_scalar( directivity, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( directivity )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "directivity_%s_fill_%.2f_att_%.1f_dB", directivity.obliquity_factor, directivity.width_over_pitch * 1e2, directivity.attenuation_dB );

        end % function str_out = string_scalar( directivity )

    end % methods (Access = protected, Hidden)

end % classdef directivity < f_numbers.f_number
