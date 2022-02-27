%
% superclass for all GPU-based platforms
%
% author: Martin F. Schiffner
% date: 2022-02-19
% modified: 2022-02-19
%
classdef gpu < platforms.platform

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)

        % independent properties
        index ( 1, 1 ) double { mustBeNonnegative, mustBeInteger, mustBeNonempty } = 0 % index of GPU

    end % properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = gpu( indices )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at most one argument
            narginchk( 0, 1 );

            % ensure existence of nonempty indices
            if nargin < 1 || isempty( indices )
                indices = 0;
            end

            % property validation functions ensure valid indices

            %--------------------------------------------------------------
            % 2.) create Wiener filter-based normalizations
            %--------------------------------------------------------------
            % constructor of superclass
            objects@platforms.platform( size( indices ) );

            % iterate normalizations
            for index_object = 1:numel( indices )

                % set independent properties
                objects( index_object ).index = indices( index_object );

            end % for index_object = 1:numel( indices )

        end % function objects = gpu( indices )

    end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( gpu )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class platforms.platform for gpu (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "gpu_%d", gpu.index );

        end % function str_out = string_scalar( gpu )

    end % methods (Access = protected, Hidden)

end % classdef gpu < platforms.platform
