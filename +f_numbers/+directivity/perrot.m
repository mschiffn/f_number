%
% superclass for all directivity-derived F-numbers
% boundary condition: soft baffle
%
% settings proposed by Perrot et al.
%
% author: Martin F. Schiffner
% date: 2021-08-06
% modified: 2022-01-08
%
classdef perrot < f_numbers.directivity.soft

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = perrot( widths_over_pitch )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % superclass ensures valid widths_over_pitch

            %--------------------------------------------------------------
            % 2.) create directivity-derived F-numbers (soft baffle)
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.directivity.soft( widths_over_pitch, 3 )

        end % function objects = perrot( widths_over_pitch )

    end % methods

end % classdef perrot < f_numbers.directivity.soft
