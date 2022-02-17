%
% abstract superclass for all types of normalization
%
% author: Martin F. Schiffner
% date: 2022-02-05
% modified: 2022-02-05
%
classdef (Abstract) normalization

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %------------------------------------------------------------------
        % constructor
        %------------------------------------------------------------------
        function objects = normalization( size )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % calling method ensures existence of nonempty size

            % ensure row vector for size
            if ~isrow( size )
                errorStruct.message = 'size must be a row vector!';
                errorStruct.identifier = 'normalization:NoRowVector';
                error( errorStruct );
            end

            % ensure nonempty positive integers
            mustBePositive( size );
            mustBeInteger( size );
            mustBeNonempty( size );

            %--------------------------------------------------------------
            % 2.) create windows
            %--------------------------------------------------------------
            % repeat default normalization
            objects = repmat( objects, size );

        end % function objects = normalization( size )

        %------------------------------------------------------------------
        % apply normalization
        %------------------------------------------------------------------
        function samples_norm = apply( normalizations, samples )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure two arguments
            narginchk( 2, 2 );

            % ensure class normalizations.normalization
            if ~isa( normalizations, 'normalizations.normalization' )
                errorStruct.message = 'normalizations must be normalizations.normalization!';
                errorStruct.identifier = 'apply:NoNormalizations';
                error( errorStruct );
            end

            % ensure cell array for samples
            if ~iscell( samples )
                samples = { samples };
            end

            % ensure equal number of dimensions and sizes
            [ normalizations, samples ] = auxiliary.ensureEqualSize( normalizations, samples );

            %--------------------------------------------------------------
            % 2.) compute samples
            %--------------------------------------------------------------
            % specify cell array for samples_norm
            samples_norm = cell( size( normalizations ) );

            % iterate normalizations
            for index_object = 1:numel( normalizations )

                % compute samples (scalar)
                samples_norm{ index_object } = apply_scalar( normalizations( index_object ), samples{ index_object } );

            end % for index_object = 1:numel( normalizations )

            % avoid cell array for scalar normalizations
            if isscalar( normalizations )
                samples_norm = samples_norm{ 1 };
            end

        end % function samples_norm = apply( normalizations, samples )

        %------------------------------------------------------------------
        % string array (overload string method)
        %------------------------------------------------------------------
        function strs_out = string( normalizations )

            %--------------------------------------------------------------
            % 1.) check arguments
            %--------------------------------------------------------------
            % ensure one argument
            narginchk( 1, 1 );

            % ensure class normalizations.normalization
            if ~isa( normalizations, 'normalizations.normalization' )
                errorStruct.message = 'normalizations must be normalizations.normalization!';
                errorStruct.identifier = 'string:NoNormalizations';
                error( errorStruct );
            end

            %--------------------------------------------------------------
            % 2.) string array
            %--------------------------------------------------------------
            % repeat empty string
            strs_out = repmat( "", size( normalizations ) );

            % iterate normalizations
            for index_object = 1:numel( normalizations )

                strs_out( index_object ) = string_scalar( normalizations( index_object ) );

            end % for index_object = 1:numel( normalizations )

        end % function strs_out = string( normalizations )

    end % methods

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% methods (Abstract, protected, hidden)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Abstract, Access = protected, Hidden)

        %------------------------------------------------------------------
        % apply normalization (scalar)
        %------------------------------------------------------------------
        samples = apply_scalar( normalization, samples );

        %------------------------------------------------------------------
        % string array (scalar)
        %------------------------------------------------------------------
        str_out = string_scalar( normalization );

    end % methods (Abstract, Access = protected, Hidden)

end % classdef (Abstract) normalization
