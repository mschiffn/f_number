%
% superclass for all Tukey (cosine-tapered) windows
%
% author: Martin F. Schiffner
% date: 2021-08-10
% modified: 2023-12-20
%
classdef tukey < windows.window

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess = private)

        % independent properties
        fraction_cosine ( 1, 1 ) double { mustBePositive, mustBeLessThan( fraction_cosine, 1 ), mustBeNonempty } = 0.5 % cosine fraction

        % dependent properties
        fraction_rectangle = 0.5 % rectangle fraction

	end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = tukey( fractions_cosine )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at most one argument
            narginchk( 0, 1 );

            % ensure existence of nonempty fractions_cosine
            if nargin < 1 || isempty( fractions_cosine )
                fractions_cosine = 0.5;
            end

            % use boxcar windows
            if any( fractions_cosine <= 0, 'all' )
                errorStruct.message = 'Use windows.boxcar for a cosine fraction of zero!';
                errorStruct.identifier = 'tukey:InvalidFractionsCosine';
                error( errorStruct );
            end

            % use Hann windows
            if any( fractions_cosine >= 1, 'all' )
                errorStruct.message = 'Use windows.hann for a cosine fraction of one!';
                errorStruct.identifier = 'tukey:InvalidFractionsCosine';
                error( errorStruct );
            end

            % property validation functions ensure valid fractions_cosine

            %--------------------------------------------------------------
            % 2.) create Tukey windows
            %--------------------------------------------------------------
            % constructor of superclass
            objects@windows.window( size( fractions_cosine ) );

            % iterate Tukey windows
            for index_object = 1:numel( fractions_cosine )

                % set independent properties
                objects( index_object ).fraction_cosine = fractions_cosine( index_object );

                % set dependent properties
                objects( index_object ).fraction_rectangle = 1 - objects( index_object ).fraction_cosine;

            end % for index_object = 1:numel( fractions_cosine )

        end % function objects = tukey( fractions_cosine )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute samples (scalar)
        %------------------------------------------------------------------
        function samples = compute_samples_scalar( tukey, positions_over_halfwidth_abs )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window for tukey (scalar)
            % calling method ensures positions_over_halfwidth_abs < 1

            %--------------------------------------------------------------
            % 2.) compute samples (scalar)
            %--------------------------------------------------------------
            % position indicators
            indicator_inside = positions_over_halfwidth_abs < 1;
            samples = double( indicator_inside );
            indicator_taper = indicator_inside & ( positions_over_halfwidth_abs > tukey.fraction_rectangle );
            positions_over_halfwidth_abs_diff_over_length = ( positions_over_halfwidth_abs( indicator_taper ) - tukey.fraction_rectangle ) ./ tukey.fraction_cosine;
            samples( indicator_taper ) = ( 1 + cos( pi * positions_over_halfwidth_abs_diff_over_length ) ) / 2;

        end % function samples = compute_samples_scalar( tukey, positions_over_halfwidth_abs )

        %------------------------------------------------------------------
        % compute derivatives (scalar)
        %------------------------------------------------------------------
        function derivatives = compute_derivatives_scalar( tukey, positions_over_halfwidth_abs )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for tukey (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute derivatives (scalar)
            %--------------------------------------------------------------
            % compute lower and upper bounds
            positions_over_halfwidth_abs_thresh = 1 - tukey.fraction_cosine;

            % absolute values of the positions
            positions_over_halfwidth_abs = abs( positions_over_halfwidth_abs );

            % value of first derivative
            derivatives = zeros( size( positions_over_halfwidth_abs ) );
            indicator_taper = ( positions_over_halfwidth_abs > positions_over_halfwidth_abs_thresh ) & ( positions_over_halfwidth_abs < 1 );
            positions_over_halfwidth_abs_diff_over_length = ( positions_over_halfwidth_abs - positions_over_halfwidth_abs_thresh ) ./ tukey.fraction_cosine;
            derivatives( indicator_taper ) = -pi * sign( positions_over_halfwidth_abs ) * sin( pi * positions_over_halfwidth_abs_diff_over_length( indicator_taper ) ) / ( 2 * tukey.fraction_cosine );

        end % function derivatives = compute_derivatives_scalar( tukey, positions_over_halfwidth_abs )

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
            str_out = sprintf( "tukey_%.2f", tukey.fraction_cosine );

        end % function str_out = string_scalar( tukey )

	end % methods (Access = protected, Hidden)

end % classdef tukey < windows.window
