%
% superclass for all CPU-based platforms
%
% author: Martin F. Schiffner
% date: 2022-02-19
% modified: 2022-02-19
%
classdef cpu < platforms.platform

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = cpu( size )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at most one argument
            narginchk( 0, 1 );

            % ensure existence of nonempty size
            if nargin < 1 || isempty( size )
                size = [ 1, 1 ];
            end

            %--------------------------------------------------------------
            % 2.) create inactive normalizations
            %--------------------------------------------------------------
            % constructor of superclass
            objects@platforms.platform( size );

		end % function objects = cpu( size )

    end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % apply normalization (scalar)
        %------------------------------------------------------------------
        function samples_norm = apply_scalar( ~, samples )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class platforms.platform for wiener (scalar)
            % calling method ensures for samples

            %--------------------------------------------------------------
            % 2.) compute samples (scalar)
            %--------------------------------------------------------------
            samples_sum = sum( samples, 2 );
            samples_norm = samples ./ samples_sum;

        end % function samples_norm = apply_scalar( ~, samples )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( ~ )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class platforms.platform for cpu (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "cpu" );

        end % function str_out = string_scalar( ~ )

    end % methods (Access = protected, Hidden)

end % classdef cpu < platforms.platform

