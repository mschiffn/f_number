%
% superclass for all directivity-derived F-numbers according to
% Szabo (hard baffle model)
%
% author: Martin F. Schiffner
% date: 2021-08-06
% modified: 2021-08-19
%
classdef szabo < f_numbers.directivity.directivity

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = szabo( varargin )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % superclass ensures valid varargin

            %--------------------------------------------------------------
            % 2.) create directivity-derived F-numbers (Szabo)
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.directivity.directivity( varargin{ : } );

        end % function objects = szabo( varargin )

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

            %
            values = element_width_over_lambda / 1.2;

%             % directivity pattern (hard baffle)
%             directivity = @( theta, width_norm ) sinc( width_norm * sin( theta ) );
% 
%             % function to minimize
%             f = @( theta, width_norm ) abs( directivity( theta, width_norm ) - szabo.thresh );
% 
%             % initialize F-numbers w/ zeros
%             values = zeros( size( element_pitch_over_lambda ) );
% 
%             % iterate samples
%             for index_width = 1:numel( element_pitch_over_lambda )
%                 alpha = fminbnd( @( theta ) f( theta, element_width_over_lambda( index_width ) ), 0, pi / 2 );
%                 values( index_width ) = 1 / ( 2 * tan( alpha ) );
%             end

        end % function values = compute_values_scalar( szabo, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( szabo )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "szabo_fill_%.2f_att_%.1f", szabo.width_over_pitch * 1e2, szabo.attenuation_dB );

        end % function str_out = string_scalar( szabo )

	end % methods (Access = protected, Hidden)

end % classdef szabo < f_numbers.directivity.directivity
