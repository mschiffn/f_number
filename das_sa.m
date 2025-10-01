function [ image, F_number_tx_values, F_number_rx_values, signal ] = das_sa( positions_x, positions_z, data_RF, f_s, element_width, element_pitch, c_0, f_bounds, index_t0, window_tx, window_rx, F_number_tx, F_number_rx, normalization_tx, normalization_rx, platform )
% DAS_SA Delay-and-Sum (DAS) Beamforming [ Fourier domain, synthetic aperture ]
%
% Form an ultrasound image from
% a complete synthetic aperture scan by using
% the delay-and-sum (DAS) algorithm in
% the Fourier domain.
%
% The Fourier domain permits:
%  1.) independent focusing of different frequencies;
%  2.) exact corrections of the arrival times independent of the sampling rate;
%  3.) elimination of out-of-band noise; and
%  4.) usage of frequency-dependent apodization weights.
%
% A uniform linear transducer array is assumed.
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% minimal:
% image = das_sa( positions_x, positions_z, data_RF, f_s, element_width, element_pitch, c_0 );
%
% maximal:
% [ image, F_number_tx_values, F_number_rx_values, signal ] = das_sa( positions_x, positions_z, data_RF, f_s, element_width, element_pitch, c_0, f_bounds, index_t0, window_tx, window_rx, F_number_tx, F_number_rx, normalization_tx, normalization_rx, platform );
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% REQUIRED
%   01.) positions_x:       lateral voxel positions (m)
%   02.) positions_z:       axial voxel positions (m)
%   03.) data_RF:           RF data of a full SA scan (3d array; 1st dimension: time, 2nd dimension: rx array element index, 3rd dimension: tx array element index)
%   04.) f_s:               sampling rate of the RF data (Hz)
%   05.) element_width:     element width of the uniform linear array (m)
%   06.) element_pitch:     element pitch of the uniform linear array (m)
%   07.) c_0:               speed of sound (m/s)
%
% OPTIONAL
%   08.) f_bounds:          frequency bounds (Hz)
%   09.) index_t0:          time index of the sample extracted from the focused RF signal (1)
%   10.) window_tx:         window functions for transmit apodization ( object of class windows.window )
%   11.) window_rx:         window functions for receive apodization ( object of class windows.window )
%   12.) F_number_tx:       transmit F-number ( object of class f_numbers.f_number )
%   13.) F_number_rx:       receive F-number ( object of class f_numbers.f_number )
%   14.) normalization_tx:  normalization of the complex-valued transmit apodization weights ( object of class normalizations.normalization )
%   15.) normalization_rx:  normalization of the complex-valued receive apodization weights ( object of class normalizations.normalization )
%   16.) platform:          platform for all computations ( object of class platforms.platform )
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%   01.) image:               complex-valued DAS image
%   02.) F_number_tx_values:  value of the transmit F-number for each frequency
%   03.) F_number_rx_values:  value of the receive F-number for each frequency
%   04.) signal:              focused RF signal for last image voxel
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%   [1] M. F. Schiffner and G. Schmitz,
%       "Frequency-dependent F-number suppresses grating lobes and improves the lateral resolution in coherent plane-wave compounding,"
%       IEEE Trans. Ultrason., Ferroelectr., Freq. Control (2023).
%       DOI: <a href="matlab:web('https://doi.org/10.1109/TUFFC.2023.3291612')">10.1109/TUFFC.2023.3291612</a>
%
%   [2] M. F. Schiffner and G. Schmitz,
%       "Frequency-dependent F-number increases the contrast and the spatial resolution in fast pulse-echo ultrasound imaging,"
%       2021 IEEE Int. Ultrasonics Symp. (IUS), Xi’an, China, Sep. 2021, pp. 1–4.
%       DOI: <a href="matlab:web('https://doi.org/10.1109/IUS52206.2021.9593488')">10.1109/IUS52206.2021.9593488</a>
%       arxiv: <a href="matlab:web('https://arxiv.org/abs/2111.04593')">2111.04593</a>
%       YouTube: <a href="matlab:web('https://www.youtube.com/watch?v=T6BoYRvQ6rg')">T6BoYRvQ6rg</a>
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% - The steering of the pulse-echo beam
%   (i.e., lateral voxel positions outside the lateral array bounds) is
%   supported but not recommended.
%
% -------------------------------------------------------------------------
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   date: 2025-04-02
%   modified: 2025-10-01

%--------------------------------------------------------------------------
% 0.) check arguments
%--------------------------------------------------------------------------
% ensure at least 7 and at most 16 arguments
narginchk( 7, 16 );

% ensure at least 1 and at most 4 output arguments
nargoutchk( 1, 4 );

% ensure numeric and real-valued row vector for positions_x
if ~( isnumeric( positions_x ) && isreal( positions_x ) && isrow( positions_x ) )
    errorStruct.message = 'positions_x must be a numeric and real-valued row vector!';
    errorStruct.identifier = 'das_sa:NoNumericAndRealRowVector';
    error( errorStruct );
end

% ensure finite x-coordinates
if ~all( isfinite( positions_x ) )
    errorStruct.message = 'positions_x must contain finite x-coordinates!';
    errorStruct.identifier = 'das_sa:NoFiniteXCoordinates';
    error( errorStruct );
end

% ensure numeric and real-valued row vector for positions_z
if ~( isnumeric( positions_z ) && isreal( positions_z ) && isrow( positions_z ) )
    errorStruct.message = 'positions_z must be a numeric and real-valued row vector!';
    errorStruct.identifier = 'das_sa:NoNumericAndRealRowVector';
    error( errorStruct );
end

% ensure positive and finite z-coordinates
if ~( all( positions_z > 0 ) && all( isfinite( positions_z ) ) )
    errorStruct.message = 'positions_z must contain positive and finite z-coordinates!';
    errorStruct.identifier = 'das_sa:NoPositiveAndFiniteZCoordinates';
    error( errorStruct );
end

% % ensure numeric and real-valued matrix for data_RF
% if ~( isnumeric( data_RF ) && isreal( data_RF ) && ismatrix( data_RF ) )
%     errorStruct.message = 'data_RF must be a numeric and real-valued matrix!';
%     errorStruct.identifier = 'das_sa:NoNumericAndRealMatrix';
%     error( errorStruct );
% end

% ensure at least two array elements and two samples
if ~( size( data_RF, 1 ) > 1 && size( data_RF, 2 ) > 1 )
    errorStruct.message = 'data_RF must contain at least two temporal and two spatial samples!';
    errorStruct.identifier = 'das_sa:NoSamples';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for f_s
if ~( isnumeric( f_s ) && isreal( f_s ) && isscalar( f_s ) )
    errorStruct.message = 'f_s must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_sa:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure positive and finite f_s
if ~( f_s > 0 && f_s < Inf )
    errorStruct.message = 'f_s must be positive and finite!';
    errorStruct.identifier = 'das_sa:InvalidSamplingRate';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for element_width
if ~( isnumeric( element_width ) && isreal( element_width ) && isscalar( element_width ) )
    errorStruct.message = 'element_width must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_sa:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure positive and finite element_width
if ~( element_width > 0 && element_width < Inf )
    errorStruct.message = 'element_width must be positive and finite!';
    errorStruct.identifier = 'das_sa:InvalidElementWidth';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for element_pitch
if ~( isnumeric( element_pitch ) && isreal( element_pitch ) && isscalar( element_pitch ) )
    errorStruct.message = 'element_pitch must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_sa:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure larger element_pitch than element_width and finite element_pitch
if ~( element_pitch > element_width && element_pitch < Inf )
    errorStruct.message = 'element_pitch must be larger than element_width and finite!';
    errorStruct.identifier = 'das_sa:InvalidElementPitch';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for c_0
if ~( isnumeric( c_0 ) && isreal( c_0 ) && isscalar( c_0 ) )
    errorStruct.message = 'c_0 must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_sa:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure positive and finite c_0
if ~( c_0 > 0 && c_0 < Inf )
    errorStruct.message = 'c_0 must be positive and finite!';
    errorStruct.identifier = 'das_sa:InvalidSpeedOfSound';
    error( errorStruct );
end

% ensure existence of nonempty f_bounds
if nargin < 8 || isempty( f_bounds )
    f_bounds = [ eps, f_s / 2 ];
end

% ensure numeric and real-valued row vector for f_bounds
if ~( isnumeric( f_bounds ) && isreal( f_bounds ) && isrow( f_bounds ) )
    errorStruct.message = 'f_bounds must be a numeric and real-valued row vector!';
    errorStruct.identifier = 'das_sa:NoNumericAndRealRowVector';
    error( errorStruct );
end

% ensure valid vector for f_bounds
if ~( numel( f_bounds ) == 2 && f_bounds( 1 ) > 0 && f_bounds( 2 ) > f_bounds( 1 ) && f_bounds( 2 ) <= f_s / 2 )
    errorStruct.message = 'f_bounds must be a strictly monotonic increasing sequence!';
    errorStruct.identifier = 'das_sa:InvalidFrequencyBounds';
    error( errorStruct );
end

% ensure existence of nonempty index_t0
if nargin < 9 || isempty( index_t0 )
	index_t0 = 0;
end

% ensure numeric and real-valued scalar for index_t0
if ~( isnumeric( index_t0 ) && isreal( index_t0 ) && isscalar( index_t0 ) )
    errorStruct.message = 'index_t0 must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_sa:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure nonnegative integer for index_t0
if ~( index_t0 >= 0 && index_t0 == floor( index_t0 ) )
    errorStruct.message = 'index_t0 must be a nonnegative integer!';
    errorStruct.identifier = 'das_sa:InvalidTimeIndex';
    error( errorStruct );
end

% ensure existence of nonempty window
if nargin < 10 || isempty( window_tx )
    window_tx = windows.tukey( 0.2 );
end

% ensure class windows.window
if ~( isa( window_tx, 'windows.window' ) && isscalar( window_tx ) )
    errorStruct.message = 'window_tx must be scalar windows.window!';
    errorStruct.identifier = 'das_sa:NoWindowTX';
    error( errorStruct );
end

% ensure existence of nonempty F_number_tx
if nargin < 12 || isempty( F_number_tx )
    F_number_tx = f_numbers.grating.angle_lb( 60, 1.5 );
end

% ensure class f_numbers.f_number
if ~( isa( F_number_tx, 'f_numbers.f_number' ) && isscalar( F_number_tx ) )
	errorStruct.message = 'F_number_tx must be scalar f_numbers.f_number!';
	errorStruct.identifier = 'das_sa:NoFNumberTX';
	error( errorStruct );
end

% ensure existence of nonempty normalization_tx
if nargin < 14 || isempty( normalization_tx )
    normalization_tx = normalizations.on;
end

% ensure class normalizations.normalization
if ~( isa( normalization_tx, 'normalizations.normalization' ) && isscalar( normalization_tx ) )
    errorStruct.message = 'normalization_tx must be scalar normalizations.normalization!';
    errorStruct.identifier = 'das_sa:NoNormalizationTX';
    error( errorStruct );
end

% ensure existence of nonempty platform
if nargin < 16 || isempty( platform )
    platform = platforms.gpu( 0 );
end

% ensure class platforms.platform
if ~( isa( platform, 'platforms.platform' ) && isscalar( platform ) )
    errorStruct.message = 'platform must be scalar platforms.platform!';
    errorStruct.identifier = 'das_sa:NoPlatform';
    error( errorStruct );
end

%--------------------------------------------------------------------------
% 1.) check platform
%--------------------------------------------------------------------------
% check platform
if isa( platform, 'platforms.gpu' )
    % execute GPU code
    [ image, weights_tx, weights_rx ] = cuda.gpu_bf_saft_rf( positions_x, positions_z, data_RF, f_s, element_width, element_pitch, c_0, f_bounds, index_t0, window_tx, window_rx, F_number_tx, F_number_rx, normalization_tx, normalization_rx, platform.index );
    return
end

%--------------------------------------------------------------------------
% 2.) geometry
%--------------------------------------------------------------------------
% print status
time_start = tic;
str_date_time = sprintf( '%04d-%02d-%02d: %02d:%02d:%02d', fix( clock ) );
fprintf( '\t %s: Delay-and-Sum (DAS) Beamforming [ Fourier domain, synthetic aperture ]... ', str_date_time );

% number of array elements
N_elements = size( data_RF, 2 );
M_elements = ( N_elements - 1 ) / 2;

% half-width of the elements
element_width_over_two = element_width / 2;

% centroids and bounds of element faces
positions_ctr_x = (-M_elements:M_elements) * element_pitch;
positions_lbs_x = positions_ctr_x.' - element_width_over_two;
positions_ubs_x = positions_ctr_x.' + element_width_over_two;

% lateral distances [ N_pos_x, N_elements ]
dist_lateral = positions_ctr_x - positions_x.';

% determine left part of the aperture
indicator_left = dist_lateral < 0;

% incident wave travel times
% TODO: which element fires first;
% t_in_plus_shift = sqrt( ( positions_x - position_src( 1 ) ).^2 + ( positions_z.' - position_src( 2 ) ).^2 ) / c_0 + index_t0 / f_s;

% scattered wave travel times
t_sc_lateral_squared = ( dist_lateral / c_0 ).^2;
t_sc_axial_squared = ( positions_z / c_0 ).^2;

% maximum relative time shift in the electronic focusing ( zero padding in DFT )
N_samples_t_pad = ceil( ( sqrt( diff( positions_ctr_x( [ 1, end ] ) )^2 + positions_z( 1 )^2 ) - positions_z( 1 ) ) * f_s / c_0 );

%--------------------------------------------------------------------------
% 3.) create frequency axis
%--------------------------------------------------------------------------
% number of samples in time domain
N_samples_t = size( data_RF, 1 );

% ensure odd number of points in DFT
N_points_dft = N_samples_t + N_samples_t_pad + 1 - mod( N_samples_t + N_samples_t_pad, 2 );

% boundary frequency indices
index_f_lb = ceil( f_bounds( 1 ) * N_points_dft / f_s ) + 1;
index_f_ub = floor( f_bounds( 2 ) * N_points_dft / f_s ) + 1;
indices_f = (index_f_lb:index_f_ub).';

% frequency axes
axis_f_bp = f_s * ( indices_f - 1 ) / N_points_dft;
axis_omega_bp = 2 * pi * axis_f_bp;

%--------------------------------------------------------------------------
% 4.) frequency-dependent F-numbers
%--------------------------------------------------------------------------
% element pitch-to-wavelength ratio
element_pitch_over_lambda = element_pitch * axis_f_bp / c_0;

% compute values of the F-number for each frequency
F_number_tx_values = compute_values( F_number_tx, element_pitch_over_lambda );
F_number_rx_values = compute_values( F_number_rx, element_pitch_over_lambda );

% % maximum F-numbers for each axial position
% F_ub = sqrt( axis_f_bp .* positions_z / ( 2.88 * c_0 ) );
% 
% % apply maximum F-numbers
% F_number_values = min( F_number_values, F_ub );

%--------------------------------------------------------------------------
% 5.) receive focusing
%--------------------------------------------------------------------------
% compute DFT of RF data
data_RF_dft = fft( data_RF, N_points_dft, 1 );
data_RF_dft_analy = 2 * data_RF_dft( indices_f, :, : );

% desired half-widths of the receive aperture
width_aperture_tx_over_two_desired = positions_z ./ ( 2 * F_number_tx_values );
width_aperture_rx_over_two_desired = positions_z ./ ( 2 * F_number_rx_values );

% initialize image w/ zeros
image = complex( zeros( numel( positions_z ), numel( positions_x ) ) );

% iterate lateral voxel positions
for index_pos_x = 1:numel( positions_x )

    % print progress in percent
    fprintf( '%5.1f %%', ( index_pos_x - 1 ) / numel( positions_x ) * 1e2 );

    %----------------------------------------------------------------------
    % transmit focusing
    %----------------------------------------------------------------------
    % map desired bounds of the transmit aperture to element indices
    indices_aperture_tx_lb = max( ceil( M_elements + ( positions_x( index_pos_x ) - width_aperture_tx_over_two_desired ) / element_pitch ) + 1, 1 );
    indices_aperture_tx_ub = min( floor( M_elements + ( positions_x( index_pos_x ) + width_aperture_tx_over_two_desired ) / element_pitch ) + 1, N_elements );

    % determine frequencies to process for each axial voxel position
    indicator_valid_tx = indices_aperture_tx_lb <= indices_aperture_tx_ub;

    % half-widths of the left transmit aperture
    width_aperture_tx_left_over_two = zeros( numel( indices_f ), numel( positions_z ) );
    width_aperture_tx_left_over_two( indicator_valid_tx ) = positions_x( index_pos_x ) - positions_lbs_x( indices_aperture_tx_lb( indicator_valid_tx ) );

    % half-widths of the right transmit aperture
    width_aperture_tx_right_over_two = zeros( numel( indices_f ), numel( positions_z ) );
    width_aperture_tx_right_over_two( indicator_valid_tx ) = positions_ubs_x( indices_aperture_tx_ub( indicator_valid_tx ) ) - positions_x( index_pos_x );

    %----------------------------------------------------------------------
    % receive focusing
    %----------------------------------------------------------------------
    % map desired bounds of the transmit aperture to element indices
    indices_aperture_rx_lb = max( ceil( M_elements + ( positions_x( index_pos_x ) - width_aperture_rx_over_two_desired ) / element_pitch ) + 1, 1 );
    indices_aperture_rx_ub = min( floor( M_elements + ( positions_x( index_pos_x ) + width_aperture_rx_over_two_desired ) / element_pitch ) + 1, N_elements );

    % determine frequencies to process for each axial voxel position
    indicator_valid_rx = indices_aperture_rx_lb <= indices_aperture_rx_ub;

    % half-widths of the left transmit aperture
    width_aperture_rx_left_over_two = zeros( numel( indices_f ), numel( positions_z ) );
    width_aperture_rx_left_over_two( indicator_valid_rx ) = positions_x( index_pos_x ) - positions_lbs_x( indices_aperture_rx_lb( indicator_valid_rx ) );

    % half-widths of the right transmit aperture
    width_aperture_rx_right_over_two = zeros( numel( indices_f ), numel( positions_z ) );
    width_aperture_rx_right_over_two( indicator_valid_rx ) = positions_ubs_x( indices_aperture_rx_ub( indicator_valid_rx ) ) - positions_x( index_pos_x );

    % determine axial voxel positions that require processing
    indices_pos_z = find( any( indicator_valid_tx & indicator_valid_rx, 1 ) );

    % iterate axial voxel positions that require processing
    for index_pos_z = indices_pos_z

        %------------------------------------------------------------------
        % a) select frequencies to process
        %------------------------------------------------------------------
        indicator_valid_tx_act = indicator_valid_tx( :, index_pos_z );
        indicator_valid_rx_act = indicator_valid_rx( :, index_pos_z );

        %------------------------------------------------------------------
        % b) compute time-of-flight (TOF)
        %------------------------------------------------------------------
        t_sc = sqrt( t_sc_lateral_squared( index_pos_x, : ) + t_sc_axial_squared( index_pos_z ) ) + 0.5 * index_t0 / f_s;
        weights = exp( 1j * axis_omega_bp( indicator_valid_rx_act ) .* t_sc );

        %------------------------------------------------------------------
        % c) compute apodization weights for current voxel
        %------------------------------------------------------------------
        % check type of window function
        if isa( window_tx, 'windows.boxcar' )

            %--------------------------------------------------------------
            % i.) simple solution for boxcar window
            %--------------------------------------------------------------
            window_tx_samples = compute_samples( window_tx, dist_lateral( index_pos_x, : ) ./ width_aperture_tx_over_two_desired( indicator_valid_tx_act, index_pos_z ) );

        else

            %--------------------------------------------------------------
            % ii.) complex solution for other windows
            %--------------------------------------------------------------
            % sample window functions
            indicator_left_act = indicator_left( index_pos_x, : );
            window_tx_samples_left = compute_samples( window_tx, dist_lateral( index_pos_x, indicator_left_act ) ./ width_aperture_tx_left_over_two( indicator_valid_tx_act, index_pos_z ) );
            window_tx_samples_right = compute_samples( window_tx, dist_lateral( index_pos_x, ~indicator_left_act ) ./ width_aperture_tx_right_over_two( indicator_valid_tx_act, index_pos_z ) );
            window_tx_samples = [ window_tx_samples_left, window_tx_samples_right ];

        end % if isa( window_tx, 'windows.boxcar' )

        % check type of window function
        if isa( window_rx, 'windows.boxcar' )

            %--------------------------------------------------------------
            % i.) simple solution for boxcar window
            %--------------------------------------------------------------
            window_rx_samples = compute_samples( window_rx, dist_lateral( index_pos_x, : ) ./ width_aperture_rx_over_two_desired( indicator_valid_rx_act, index_pos_z ) );

        else

            %--------------------------------------------------------------
            % ii.) complex solution for other windows
            %--------------------------------------------------------------
            % sample window functions
            window_rx_samples_left = compute_samples( window_rx, dist_lateral( index_pos_x, indicator_left_act ) ./ width_aperture_rx_left_over_two( indicator_valid_rx_act, index_pos_z ) );
            window_rx_samples_right = compute_samples( window_rx, dist_lateral( index_pos_x, ~indicator_left_act ) ./ width_aperture_rx_right_over_two( indicator_valid_rx_act, index_pos_z ) );
            window_rx_samples = [ window_rx_samples_left, window_rx_samples_right ];

        end % if isa( window_tx, 'windows.boxcar' )

        % normalize window samples
        window_tx_samples = apply( normalization_tx, window_tx_samples );
        window_rx_samples = apply( normalization_rx, window_rx_samples );

        %------------------------------------------------------------------
        % d) apply weights and focus RF signals
        %------------------------------------------------------------------
        data_RF_dft_analy_focused_tx = sum( reshape( window_tx_samples .* weights, [ size( weights, 1 ), 1, N_elements ] ) .* data_RF_dft_analy( indicator_valid_tx_act, :, : ), 3 );
        data_RF_dft_analy_focused = sum( window_rx_samples .* weights .* data_RF_dft_analy_focused_tx, 2 );

        %------------------------------------------------------------------
        % e) inverse DFT of focused signal at time index 0
        %------------------------------------------------------------------
        image( index_pos_z, index_pos_x ) = sum( data_RF_dft_analy_focused, 1 );

    end % for index_pos_z = indices_pos_z

	% erase progress in percent
	fprintf( '\b\b\b\b\b\b\b' );

end % for index_pos_x = 1:numel( positions_x )

% return focused RF signal
if nargout >= 4
	signal = zeros( N_points_dft, 1 );
	signal( indices_f( indicator_valid_tx_act ) ) = data_RF_dft_analy_focused;
    signal = ifft( signal, N_points_dft, 1 );
end

% infer and print elapsed time
time_elapsed = toc( time_start );
fprintf( 'done! (%f s)\n', time_elapsed );

end % function [ image, F_number_tx_values, F_number_rx_values, signal ] = das_sa( positions_x, positions_z, data_RF, f_s, element_width, element_pitch, c_0, f_bounds, index_t0, window_tx, window_rx, F_number_tx, F_number_rx, normalization_tx, normalization_rx, platform )
