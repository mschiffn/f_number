%
% superclass for all directivity-derived F-numbers
% boundary condition: rigid baffle
%
% approximation and settings proposed by Szabo [1]
%
% REFERENCES:
%   [1] T. L. Szabo, Diagnostic Ultrasound Imaging: Inside Out, 2nd. Elsevier Academic Press, Dec. 2013.
%
% ABOUT:
%   author: Martin F. Schiffner
%   date: 2021-08-06
%   modified: 2023-04-15
%
classdef szabo < f_numbers.directivity.rigid

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = szabo( widths_over_pitch )
            % szabo Create directivity-derived F-numbers according to Szabo.
            %
            % This F-number is described in Eq. (10.11) in [1].
            %
            % objects = f_numbers.directivity.rigid.szabo( widths_over_pitch );
            %
            % The method is vectorized.
            % The size of objects equals the
            %
            % INPUTS:
            %   widths_over_pitch: array of element width over element pitch
            %
            % OUTPUTS:
            %   objects: array of f_numbers.directivity.rigid.szabo
            %
            % REFERENCES:
            %   [1] T. L. Szabo, Diagnostic Ultrasound Imaging: Inside Out, 2nd. Elsevier Academic Press, Dec. 2013.
            %

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % superclass ensures valid widths_over_pitch

            %--------------------------------------------------------------
            % 2.) create directivity-derived F-numbers (rigid baffle)
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.directivity.rigid( widths_over_pitch, 6 );

        end % function objects = szabo( widths_over_pitch )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( szabo, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            % normalized element width
            element_width_over_lambda = szabo.width_over_pitch * element_pitch_over_lambda;

            % compute values
            values = element_width_over_lambda / 1.2;

        end % function values = compute_values_scalar( szabo, element_pitch_over_lambda )

	end % methods (Access = protected, Hidden)

end % classdef szabo < f_numbers.directivity.rigid
