%
% abstract superclass for all directivity-derived F-numbers
%
% author: Martin F. Schiffner
% date: 2021-08-06
% modified: 2021-08-06
%
classdef (Abstract) directivity < f_numbers.f_number

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess = private)

        % independent properties
        width_over_pitch ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 0.918	% element width over element pitch (filling factor)
        attenuation_dB ( 1, 1 ) double { mustBePositive, mustBeNonempty } = 3           % tolerable attenuation in the directivity

        % dependent properties
        thresh ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 10^( -3 / 20 )	% threshold resulting from tolerable attenuation

	end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = directivity( widths_over_pitch, attenuations_dB )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure at least one and at most two arguments
            narginchk( 1, 2 );

            % property validation functions ensure valid widths_over_pitch

            % ensure existence of nonempty attenuations_dB
            if nargin < 2 || isempty( attenuations_dB )
                attenuations_dB = 3;
            end

            %--------------------------------------------------------------
            % 2.) create directivity-derived F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.f_number( size( widths_over_pitch ) );

            % iterate directivity-derived F-numbers
            for index_object = 1:numel( objects )

                % set independent properties
                objects( index_object ).width_over_pitch = widths_over_pitch( index_object );
                objects( index_object ).attenuation_dB = attenuations_dB( index_object );

                % set dependent properties
                objects( index_object ).thresh = 10^( -objects( index_object ).attenuation_dB / 20 );

            end % for index_object = 1:numel( objects )

        end % function objects = directivity( widths_over_pitch, attenuations_dB )

	end % methods

end % classdef (Abstract) directivity < f_numbers.f_number