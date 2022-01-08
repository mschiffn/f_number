%
% superclass for all directivity-derived F-numbers
% boundary condition: soft baffle
%
% author: Martin F. Schiffner
% date: 2021-08-06
% modified: 2022-01-08
%
classdef soft < f_numbers.directivity.directivity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = soft( widths_over_pitch, attenuations_dB )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at least one and at most two arguments
            narginchk( 1, 2 );

            % superclass ensures valid widths_over_pitch

            % ensure existence of attenuations_dB
            if nargin < 2
                attenuations_dB = [];
            end

            % superclass ensures valid attenuations_dB

            % superclass ensures equal number of dimensions and sizes

            %--------------------------------------------------------------
            % 2.) create directivity-derived F-numbers (soft baffle)
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.directivity.directivity( widths_over_pitch, f_numbers.directivity.obliquity.soft, attenuations_dB );

        end % function objects = soft( widths_over_pitch, attenuations_dB )

    end % methods

end % classdef soft < f_numbers.directivity.directivity
