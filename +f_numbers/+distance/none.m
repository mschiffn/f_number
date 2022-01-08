%
% superclass for all F-numbers that prevent
% the first-order grating lobes from overlapping with
% the main lobe
%
% author: Martin F. Schiffner
% date: 2022-01-07
% modified: 2022-01-07
%
classdef none < f_numbers.distance.constant

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = none( size )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at most one argument
            narginchk( 0, 1 );

            % ensure existence of nonempty size
            if nargin < 1 || isempty( size )
                size = [ 1, 1 ];
            end

            % superclass ensures valid size

            %--------------------------------------------------------------
            % 2.) create sampling theorem-derived F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.distance.constant( zeros( size ), inf( size ) );

        end % function objects = none( size )

    end % methods

end % classdef none < f_numbers.distance.constant
