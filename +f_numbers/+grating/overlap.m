%
% superclass for all grating lobe-derived F-numbers with
% overlap prevention
%
% author: Martin F. Schiffner
% date: 2021-08-07
% modified: 2021-08-07
%
classdef overlap < f_numbers.grating.grating

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = overlap( varargin )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % property validation functions ensure valid angles_lb

            %--------------------------------------------------------------
            % 2.) create grating lobe-derived F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.grating.grating( varargin{ : } );

        end % function objects = overlap( varargin )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( overlap, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            values = sqrt( max( element_pitch_over_lambda.^2, 0.25 ) - 0.25 );

        end % function values = compute_values_scalar( overlap, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( overlap )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "angle_lb_%.2f", overlap.angle_lb );

        end % function strs_out = string( overlap )

	end % methods (Access = protected, Hidden)

end % classdef overlap < f_numbers.grating.grating
