%
% superclass for all obliquity factors for a soft baffle
%
% author: Martin F. Schiffner
% date: 2022-01-08
% modified: 2022-01-08
%
classdef soft < f_numbers.directivity.obliquity.obliquity

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = soft( size )

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
            % 2.) create obliquity factors for a soft baffle
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.directivity.obliquity.obliquity( size );

        end % function objects = soft( size )

	end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods (protected, hidden)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( ~, theta )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.directivity.obliquity.obliquity for soft (scalar)
            % calling method ensures valid theta

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            values = cos( theta );

        end % function values = compute_values_scalar( ~, theta )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( ~ )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.directivity.obliquity.obliquity (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = "soft";

        end % function str_out = string_scalar( ~ )

    end % methods (Access = protected, Hidden)

end % classdef soft < f_numbers.directivity.obliquity.obliquity
