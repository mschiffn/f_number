%
% superclass for all constant F-numbers
%
% author: Martin F. Schiffner
% date: 2021-08-03
% modified: 2021-08-19
%
classdef constant < f_numbers.f_number

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess = private)

        % independent properties
        value ( 1, 1 ) double { mustBeNonnegative, mustBeNonempty } = 1	% value of the F-number

	end % properties

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = constant( vals )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % property validation function ensures valid vals

            %--------------------------------------------------------------
            % 2.) create constant F-numbers
            %--------------------------------------------------------------
            % constructor of superclass
            objects@f_numbers.f_number( size( vals ) );

            % iterate constant F-numbers
            for index_object = 1:numel( objects )

                % set independent properties
                objects( index_object ).value = vals( index_object );

            end % for index_object = 1:numel( objects )

        end % function objects = constant( vals )

	end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access = protected, Hidden)

        %------------------------------------------------------------------
        % compute values (scalar)
        %------------------------------------------------------------------
        function values = compute_values_scalar( f_number, element_pitch_over_lambda )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)
            % calling method ensures nonempty positive element_pitch_over_lambda for element_pitch_over_lambda

            %--------------------------------------------------------------
            % 2.) compute values (scalar)
            %--------------------------------------------------------------
            values = repmat( f_number.value, size( element_pitch_over_lambda ) );

        end % function values = compute_values_scalar( f_number, element_pitch_over_lambda )

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        function str_out = string_scalar( f_number )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures class f_numbers.f_number for f_number (scalar)

            %--------------------------------------------------------------
            % 2.) create string scalar
            %--------------------------------------------------------------
            str_out = sprintf( "constant_%.2f", f_number.value );

        end % function str_out = string_scalar( f_number )

	end % methods (Access = protected, Hidden)

end % classdef constant < f_numbers.f_number
