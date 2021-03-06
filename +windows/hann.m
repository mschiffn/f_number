%
% superclass for all Hann windows
%
% author: Martin F. Schiffner
% date: 2021-08-10
% modified: 2022-01-09
%
classdef hann < windows.window

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = hann( size )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at most one input
            narginchk( 0, 1 );

            % ensure existence of nonempty size
            if nargin < 1 || isempty( size )
                size = [ 1, 1 ];
            end

            % superclass ensures valid size

            %--------------------------------------------------------------
            % 2.) create hann apodizations
            %--------------------------------------------------------------
            % constructor of superclass
            objects@windows.window( size );

        end % function objects = hann( size )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute samples (scalar)
        %------------------------------------------------------------------
        function samples = compute_samples_scalar( ~, positions_over_halfwidth )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window for hann (scalar)
            % calling method ensures for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute samples (scalar)
            %--------------------------------------------------------------
            samples = ( abs( positions_over_halfwidth ) < 1 ) .* ( 1 + cos( pi * positions_over_halfwidth ) ) / 2;

        end % function samples = compute_samples_scalar( ~, positions_over_halfwidth )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( ~ )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class windows.window for hann (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = "hann";

        end % function str_out = string_scalar( ~ )

	end % methods (Access = protected, Hidden)

end % classdef hann < windows.window
