%
% superclass for all directivity-derived F-numbers according to
% Perrot et al.
%
% author: Martin F. Schiffner
% date: 2021-08-06
% modified: 2021-08-06
%
classdef perrot < f_numbers.directivity.directivity

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = perrot( varargin )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % superclass ensures valid widths_over_pitch

            %--------------------------------------------------------------
            % 2.) create directivity-derived F-numbers (Perrot et al.)
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.directivity.directivity( varargin{ : } );

        end % function objects = perrot( widths_over_pitch )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( perrot, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            % normalized element width
            element_width_over_lambda = perrot.width_over_pitch * element_pitch_over_lambda;

            % directivity pattern (soft baffle)
            directivity = @( theta, width_norm ) cos( theta ) * sinc( width_norm * sin( theta ) );

            % function to minimize
            f = @( theta, width_norm ) abs( directivity( theta, width_norm ) - perrot.thresh );

            % initialize F-numbers w/ zeros
            values = zeros( size( element_pitch_over_lambda ) );

            % iterate samples
            for index_width = 1:numel( element_pitch_over_lambda )
                alpha = fminbnd( @( theta ) f( theta, element_width_over_lambda( index_width ) ), 0, pi / 2 );
                values( index_width ) = 1 / ( 2 * tan( alpha ) );
            end

        end % function values = compute_values_scalar( perrot, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( perrot )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "perrot_fill_%.2f_att_%.1f", perrot.width_over_pitch * 1e2, perrot.attenuation_dB );

        end % function strs_out = string( f_numbers )

	end % methods (Access = protected, Hidden)

end % classdef perrot < f_numbers.directivity.directivity
