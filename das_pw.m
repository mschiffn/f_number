function [ image, F_number_values, signal ] = das_pw( positions_x, positions_z, data_RF, f_s, theta_incident, element_width, element_pitch, c_0, f_bounds, index_t0, window, F_number, normalization, platform )
% DAS_PW Delay-and-Sum (DAS) Beamforming [ Fourier domain, steered plane wave ]
%
% Form ultrasound images by using
% the delay-and-sum (DAS) algorithm in
% the Fourier domain.
%
% The Fourier domain permits
%  1.) independent receive focusing of different frequencies,
%  2.) exact corrections of the arrival times independent of the sampling rate, and
%  3.) usage of frequency-dependent apodization weights.
%
% A uniform linear transducer array is assumed.
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% minimal:
% image = das_pw( positions_x, positions_z, data_RF, f_s, theta_incident, element_width, element_pitch, c_0 );
%
% maximal:
% [ image, F_number_values, signal ] = das_pw( positions_x, positions_z, data_RF, f_s, theta_incident, element_width, element_pitch, c_0, f_bounds, index_t0, window, F_number, normalization, platform );
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% REQUIRED
%   01.) positions_x:       lateral voxel positions (m)
%   02.) positions_z:       axial voxel positions (m)
%   03.) data_RF:           RF data (2d array; 1st dimension: time, 2nd dimension: array element index)
%   04.) f_s:               sampling rate of the RF data (Hz)
%   05.) theta_incident:    steering angle of the incident plane wave (rad) [ ← = pi/2, ↓ = 0, → = -pi/2 ]
%   06.) element_width:     element width of the uniform linear array (m)
%   07.) element_pitch:     element pitch of the uniform linear array (m)
%   08.) c_0:               speed of sound (m/s)
%
% OPTIONAL
%   09.) f_bounds:          frequency bounds (Hz)
%   10.) index_t0:          time index of the sample extracted from the focused RF signal (1)
%   11.) window:            window function for receive apodization ( object of class windows.window )
%   12.) F_number:          receive F-number (  object of class f_numbers.f_number )
%   13.) normalization:     normalization of the complex-valued apodization weights ( object of class normalizations.normalization )
%   14.) platform:          platform for all computations ( object of class platforms.platform )
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%   01.) image:             complex-valued DAS image
%   02.) F_number_values:   value of the F-number for each frequency
%   03.) signal:            focused RF signal for last image voxel
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
% - The steering of the receive beam should be avoided.
%
% -------------------------------------------------------------------------
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   date: 2021-04-17
%   modified: 2023-07-31

%--------------------------------------------------------------------------
% 0.) check arguments
%--------------------------------------------------------------------------
% ensure at least 8 and at most 14 arguments
narginchk( 8, 14 );

% ensure numeric and real-valued row vector for positions_x
if ~( isnumeric( positions_x ) && isreal( positions_x ) && isrow( positions_x ) )
    errorStruct.message = 'positions_x must be a numeric and real-valued row vector!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealRowVector';
    error( errorStruct );
end

% ensure numeric and real-valued row vector for positions_z
if ~( isnumeric( positions_z ) && isreal( positions_z ) && isrow( positions_z ) )
    errorStruct.message = 'positions_z must be a numeric and real-valued row vector!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealRowVector';
    error( errorStruct );
end

% ensure positive z-coordinates
if ~all( positions_z > 0 )
    errorStruct.message = 'positions_z must contain positive z-coordinates!';
    errorStruct.identifier = 'das_pw:NoPositiveZCoordinates';
    error( errorStruct );
end

% ensure numeric and real-valued matrix for data_RF
if ~( ismatrix( data_RF ) && isnumeric( data_RF ) && isreal( data_RF ) )
    errorStruct.message = 'data_RF must be a numeric and real-valued matrix!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealMatrix';
    error( errorStruct );
end

% ensure at least two array elements and two samples
if ~( size( data_RF, 1 ) > 1 && size( data_RF, 2 ) > 1 )
    errorStruct.message = 'data_RF must contain at least two temporal and two spatial samples!';
    errorStruct.identifier = 'das_pw:NoSamples';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for f_s
if ~( isscalar( f_s ) && isnumeric( f_s ) && isreal( f_s ) )
    errorStruct.message = 'f_s must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure positive and finite f_s
if ~( f_s > 0 && f_s < Inf )
    errorStruct.message = 'f_s must be positive and finite!';
    errorStruct.identifier = 'das_pw:InvalidSamplingRate';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for theta_incident
if ~( isscalar( theta_incident ) && isnumeric( theta_incident ) && isreal( theta_incident ) )
    errorStruct.message = 'theta_incident must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure steering angle between -90° and 90°
if ~( theta_incident > -pi/2 && theta_incident < pi/2 )
    errorStruct.message = 'theta_incident must be larger than negative pi over two and smaller than pi over two!';
    errorStruct.identifier = 'das_pw:InvalidSteeringAngle';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for element_width
if ~( isscalar( element_width ) && isnumeric( element_width ) && isreal( element_width ) )
    errorStruct.message = 'element_width must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure positive and finite element_width
if ~( element_width > 0 && element_width < Inf )
    errorStruct.message = 'element_width must be positive and finite!';
    errorStruct.identifier = 'das_pw:InvalidElementWidth';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for element_pitch
if ~( isscalar( element_pitch ) && isnumeric( element_pitch ) && isreal( element_pitch ) )
    errorStruct.message = 'element_pitch must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure larger element_pitch than element_width and finite element_pitch
if ~( element_pitch > element_width && element_pitch < Inf )
    errorStruct.message = 'element_pitch must be larger than element_width and finite!';
    errorStruct.identifier = 'das_pw:InvalidElementPitch';
    error( errorStruct );
end

% ensure numeric and real-valued scalar for c_0
if ~( isscalar( c_0 ) && isnumeric( c_0 ) && isreal( c_0 ) )
    errorStruct.message = 'c_0 must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure positive and finite c_0
if ~( c_0 > 0 && c_0 < Inf )
    errorStruct.message = 'c_0 must be positive and finite!';
    errorStruct.identifier = 'das_pw:InvalidSpeedOfSound';
    error( errorStruct );
end

% ensure existence of nonempty f_bounds
if nargin < 9 || isempty( f_bounds )
    f_bounds = [ eps, f_s / 2 ];
end

% ensure numeric and real-valued row vector for f_bounds
if ~( isrow( f_bounds ) && isnumeric( f_bounds ) && isreal( f_bounds ) )
    errorStruct.message = 'f_bounds must be a numeric and real-valued row vector!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealRowVector';
    error( errorStruct );
end

% ensure valid vector for f_bounds
if ~( numel( f_bounds ) == 2 && f_bounds( 1 ) > 0 && f_bounds( 2 ) > f_bounds( 1 ) && f_bounds( 2 ) <= f_s / 2 )
    errorStruct.message = 'f_bounds must be a strictly monotonic increasing sequence!';
    errorStruct.identifier = 'das_pw:InvalidFrequencyBounds';
    error( errorStruct );
end

% ensure existence of nonempty index_t0
if nargin < 10 || isempty( index_t0 )
	index_t0 = 0;
end

% ensure numeric and real-valued scalar for index_t0
if ~( isscalar( index_t0 ) && isnumeric( index_t0 ) && isreal( index_t0 ) )
    errorStruct.message = 'index_t0 must be a numeric and real-valued scalar!';
    errorStruct.identifier = 'das_pw:NoNumericAndRealScalar';
    error( errorStruct );
end

% ensure nonnegative integer for index_t0
if ~( index_t0 >= 0 && index_t0 == floor( index_t0 ) )
    errorStruct.message = 'index_t0 must be a nonnegative integer!';
    errorStruct.identifier = 'das_pw:InvalidTimeIndex';
    error( errorStruct );
end

% ensure existence of nonempty window
if nargin < 11 || isempty( window )
    window = windows.tukey( 0.2 );
end

% ensure class windows.window
if ~( isa( window, 'windows.window' ) && isscalar( window ) )
    errorStruct.message = 'window must be scalar windows.window!';
    errorStruct.identifier = 'das_pw:NoWindow';
    error( errorStruct );
end

% ensure existence of nonempty F_number
if nargin < 12 || isempty( F_number )
    F_number = f_numbers.grating.angle_lb( 60, 1.5 );
end

% ensure class f_numbers.f_number
if ~( isa( F_number, 'f_numbers.f_number' ) && isscalar( F_number ) )
	errorStruct.message = 'F_number must be scalar f_numbers.f_number!';
	errorStruct.identifier = 'das_pw:NoFNumber';
	error( errorStruct );
end

% ensure existence of nonempty normalization
if nargin < 13 || isempty( normalization )
    normalization = normalizations.on;
end

% ensure class normalizations.normalization
if ~( isa( normalization, 'normalizations.normalization' ) && isscalar( normalization ) )
    errorStruct.message = 'normalization must be scalar normalizations.normalization!';
    errorStruct.identifier = 'das_pw:NoNormalization';
    error( errorStruct );
end

% ensure existence of nonempty platform
if nargin < 14 || isempty( platform )
    platform = platforms.cpu;
end

% ensure class platforms.platform
if ~( isa( platform, 'platforms.platform' ) && isscalar( platform ) )
    errorStruct.message = 'platform must be scalar platforms.platform!';
    errorStruct.identifier = 'das_pw:NoPlatform';
    error( errorStruct );
end

%--------------------------------------------------------------------------
% 1.) check platform
%--------------------------------------------------------------------------
% check platform
if isa( platform, 'platforms.gpu' )
    % execute GPU code string( normalization.window )
    % TODO: fix redefinition of theta_incident
    image = cuda.gpu_bf_das_pw_rf( positions_x, positions_z, data_RF, f_s, theta_incident + pi/2, element_width, element_pitch, f_bounds( 1 ), f_bounds( 2 ), c_0, index_t0, window, F_number, normalization, platform.index );
    return
end

%--------------------------------------------------------------------------
% 2.) geometry
%--------------------------------------------------------------------------
% print status
time_start = tic;
str_date_time = sprintf( '%04d-%02d-%02d: %02d:%02d:%02d', fix( clock ) );
fprintf( '\t %s: Delay-and-Sum (DAS) Beamforming [ Fourier domain, steered plane wave ]... ', str_date_time );

% number of array elements
N_elements = size( data_RF, 2 );
M_elements = ( N_elements - 1 ) / 2;

% half-width of the elements
element_width_over_two = element_width / 2;

% centroids and bounds of element faces
positions_ctr_x = (-M_elements:M_elements) * element_pitch;
positions_lbs_x = positions_ctr_x.' - element_width_over_two;
positions_ubs_x = positions_ctr_x.' + element_width_over_two;

% propagation direction ( ← = pi/2, ↓ = 0, → = -pi/2 )
e_theta_x = -sin( theta_incident );
e_theta_z = cos( theta_incident );

% reference position
position_ctr_x_ref = sign( e_theta_x ) * positions_ctr_x( 1 );

% lateral distances [ N_pos_x, N_elements ]
dist_lateral = positions_ctr_x - positions_x.';

% determine left part of the aperture
indicator_left = dist_lateral < 0;

% incident wave travel times
t_in_plus_shift = ( e_theta_x * ( positions_x - position_ctr_x_ref ) + e_theta_z * positions_z.' ) / c_0 + index_t0 / f_s;

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
% 4.) frequency-dependent F-number
%--------------------------------------------------------------------------
% element pitch-to-wavelength ratio
element_pitch_over_lambda = element_pitch * axis_f_bp / c_0;

% compute values of the F-number for each frequency
F_number_values = compute_values( F_number, element_pitch_over_lambda );

% maximum F-numbers for each axial position
F_ub = sqrt( axis_f_bp .* positions_z / ( 2.88 * c_0 ) );

% apply maximum F-numbers
F_number_values = min( F_number_values, F_ub );

%--------------------------------------------------------------------------
% 5.) receive focusing
%--------------------------------------------------------------------------
% compute DFT of RF data
data_RF_dft = fft( data_RF, N_points_dft, 1 );
data_RF_dft_analy = 2 * data_RF_dft( indices_f, : );

% desired half-widths of the receive subapertures
width_aperture_over_two_desired = positions_z ./ ( 2 * F_number_values );

% initialize image w/ zeros
image = zeros( numel( positions_z ), numel( positions_x ) );

% iterate lateral voxel positions
for index_pos_x = 1:numel( positions_x )

	% print progress in percent
    fprintf( '%5.1f %%', ( index_pos_x - 1 ) / numel( positions_x ) * 1e2 );

    % map desired bounds of the receive aperture to element indices
    indices_aperture_lb = max( ceil( M_elements + ( positions_x( index_pos_x ) - width_aperture_over_two_desired ) / element_pitch ) + 1, 1 );
    indices_aperture_ub = min( floor( M_elements + ( positions_x( index_pos_x ) + width_aperture_over_two_desired ) / element_pitch ) + 1, N_elements );

    % frequency-dependent processing of each z-coordinate
    indicator_valid = indices_aperture_lb <= indices_aperture_ub;

    % iterate axial voxel positions
    for index_pos_z = 1:numel( positions_z )

        %------------------------------------------------------------------
        % a) select frequencies to process
        %------------------------------------------------------------------
        indicator_valid_act = indicator_valid( :, index_pos_z );

        % actual width of the receive aperture
        width_aperture_left_over_two = positions_x( index_pos_x ) - positions_lbs_x( indices_aperture_lb( indicator_valid_act, index_pos_z ) );
        width_aperture_right_over_two = positions_ubs_x( indices_aperture_ub( indicator_valid_act, index_pos_z ) ) - positions_x( index_pos_x );

        %------------------------------------------------------------------
        % b) compute time-of-flight (TOF)
        %------------------------------------------------------------------
        t_sc = sqrt( t_sc_lateral_squared( index_pos_x, : ) + t_sc_axial_squared( index_pos_z ) );
        tof = t_in_plus_shift( index_pos_z, index_pos_x ) + t_sc;

        weights = exp( 1j * axis_omega_bp( indicator_valid_act ) * tof );

        %------------------------------------------------------------------
        % c) compute apodization weights for current voxel
        %------------------------------------------------------------------
        % check type of window function
        if isa( window, 'windows.boxcar' )

            %--------------------------------------------------------------
            % i.) simple solution for boxcar window
            %--------------------------------------------------------------
            window_samples = compute_samples( window, dist_lateral( index_pos_x, : ) ./ width_aperture_over_two_desired( indicator_valid_act, index_pos_z ) );

        else

            %--------------------------------------------------------------
            % ii.) complex solution for other windows
            %--------------------------------------------------------------
            % sample window functions
            indicator_left_act = indicator_left( index_pos_x, : );
            window_samples_left = compute_samples( window, dist_lateral( index_pos_x, indicator_left_act ) ./ width_aperture_left_over_two );
            window_samples_right = compute_samples( window, dist_lateral( index_pos_x, ~indicator_left_act ) ./ width_aperture_right_over_two );
            window_samples = [ window_samples_left, window_samples_right ];

        end % if isa( window, 'windows.boxcar' )

        % normalize window samples
        window_samples = apply( normalization, window_samples );

        %------------------------------------------------------------------
        % d) apply weights and focus RF signals
        %------------------------------------------------------------------
        data_RF_dft_analy_focused = sum( window_samples .* weights .* data_RF_dft_analy( indicator_valid_act, : ), 2 );

%         figure( index_pos_z + 3 );
%         temp = zeros( index_f_ub, N_elements );
%         temp( indices_f, : ) = window_samples .* weights;
%         window_samples_td = ifft( temp, N_points_dft, 1, 'symmetric' );
%         imagesc( illustration.dB( window_samples_td, 20 ), [ -60, 0 ] );
%         temp = zeros( index_f_ub, 1 );
%         temp( indices_f ) = data_RF_dft_analy_focused;
%         temp = ifft( temp, N_points_dft, 1, 'symmetric' );
%         plot( (1:N_samples_t), temp, (1:N_samples_t), abs( hilbert( temp ) ) );

        %------------------------------------------------------------------
        % e) inverse DFT of focused signal at time index 0
        %------------------------------------------------------------------
        image( index_pos_z, index_pos_x ) = sum( data_RF_dft_analy_focused, 1 );

    end % for index_pos_z = 1:numel( positions_z )

	% erase progress in percent
	fprintf( '\b\b\b\b\b\b\b' );

end % for index_pos_x = 1:numel( positions_x )

% return focused RF signal
if nargout >= 3
	signal = zeros( N_points_dft, 1 );
	signal( indices_f ) = data_RF_dft_analy_focused;
    signal = ifft( signal, N_points_dft, 1 );
end

% infer and print elapsed time
time_elapsed = toc( time_start );
fprintf( 'done! (%f s)\n', time_elapsed );

end % function [ image, F_number_values, signal ] = das_pw( positions_x, positions_z, data_RF, f_s, theta_incident, element_width, element_pitch, c_0, f_bounds, index_t0, window, F_number, normalization, platform )
