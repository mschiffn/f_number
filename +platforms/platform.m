%
% abstract superclass for all types of platform
%
% author: Martin F. Schiffner
% date: 2022-02-19
% modified: 2022-02-19
%
classdef (Abstract) platform

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = platform( size )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures existence of nonempty size

            % ensure row vector for size
            if ~isrow( size )
                errorStruct.message = 'size must be a row vector!';
                errorStruct.identifier = 'platform:NoRowVector';
                error( errorStruct );
            end

            % ensure nonempty positive integers
            mustBePositive( size );
            mustBeInteger( size );
            mustBeNonempty( size );

            %--------------------------------------------------------------
            % 2.) create windows
            %--------------------------------------------------------------
            % repeat default platform
            objects = repmat( objects, size );

        end % function objects = platform( size )

        %------------------------------------------------------------------
        % string array (overload string method)
        %------------------------------------------------------------------
        function strs_out = string( platforms )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % ensure class platforms.platform
            if ~isa( platforms, 'platforms.platform' )
                errorStruct.message = 'platforms must be platforms.platform!';
                errorStruct.identifier = 'string:Noplatforms';
                error( errorStruct );
            end

            %--------------------------------------------------------------
            % 2.) string array
            %--------------------------------------------------------------
            % repeat empty string
            strs_out = repmat( "", size( platforms ) );

            % iterate platforms
            for index_object = 1:numel( platforms )

                strs_out( index_object ) = string_scalar( platforms( index_object ) );

            end % for index_object = 1:numel( platforms )

        end % function strs_out = string( platforms )

    end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (Abstract, protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Abstract, Access = protected, Hidden)

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        str_out = string_scalar( normalization );

    end % methods (Abstract, Access = protected, Hidden)

end % classdef (Abstract) normalization
