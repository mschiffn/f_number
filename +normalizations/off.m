%
% superclass for all inactive normalizations
%
% author: Martin F. Schiffner
% date: 2022-02-05
% modified: 2022-02-05
%
classdef off < normalizations.normalization

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = off( size )

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
            objects@normalizations.normalization( size );

		end % function objects = off( size )

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods (protected, hidden)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % apply normalization (scalar)
        %------------------------------------------------------------------
        function samples = apply_scalar( ~, samples )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class normalizations.normalization for wiener (scalar)
            % calling method ensures for samples

            %--------------------------------------------------------------
            % 2.) apply normalization (scalar)
            %--------------------------------------------------------------
            % copy samples

        end % function samples = apply_scalar( ~, samples )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( off )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class normalizations.normalization for off (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "off" );

        end % function str_out = string_scalar( off )

    end % methods (Access = protected, Hidden)

end % classdef off < normalizations.normalization

