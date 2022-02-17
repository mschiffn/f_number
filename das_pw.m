function [ image, F_number_values, signal ] = das_pw( positions_x, positions_z, data_RF, f_s, e_theta, element_width, element_pitch, f_bounds, c_0, index_t0, window, F_number, normalization )
% DAS_PW Delay-and-Sum (DAS) Beamforming [ Fourier domain, steered plane wave ]
%
% Computes a sonogram using the delay-and-sum (DAS) algorithm in
% the Fourier domain.
%
% The Fourier domain permits
%  1.) the independent receive focusing of different frequencies,
%  2.) exact corrections of the round-trip times-of-flight independent of the sampling rate, and
%  3.) the usage of frequency-dependent apodization weights.
%
% A uniform linear transducer array is assumed.
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
%   01.) positions_x:         lateral voxel positions (m)
%   02.) positions_z:         axial voxel positions (m)
%   03.) data_RF:             RF data (2d array; 1st dimension: time, 2nd dimension: array element index)
%   04.) f_s:                 sampling rate of the RF data (Hz)
%   05.) e_theta:             propagation direction of the incident plane wave (1)
%   06.) element_width:       element width of the linear transducer array (m)
%   07.) element_pitch:       element pitch of the linear transducer array (m)
%   08.) f_bounds:            frequency bounds (Hz)
%   09.) c_0:                 speed of sound (m/s)
%   10.) index_t0:            time index of the sample extracted from the focused RF signal (1)
%   11.) window:              window function for receive apodization (  object of class windows.window )
%   12.) F_number:            receive F-number (  object of class f_numbers.f_number )
%   13.) normalization:       normalization of the complex-valued apodization weights ( object of class normalizations.normalization )
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%   01.) image:               complex-valued DAS image
%   02.) F_number_values:     value of the F-number for each frequency
%   03.) signal:              focused RF signal for last image voxel
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%   [1] M. F. Schiffner and G. Schmitz, "Frequency-Dependent F-Number Improves the Contrast and the Lateral Resolution in Coherent Plane-Wave Compounding," in press
%   [2] M. F. Schiffner and G. Schmitz, "Frequency-Dependent F-Number Increases the Contrast and the Spatial Resolution in Fast Pulse-Echo Ultrasound Imaging," in
%       2021 IEEE Int. Ultrasonics Symp. (IUS), Virtual Symposium, Sep. 2021, in press.
%       DOI:
%       arXiv: https://arxiv.org/abs/2111.04593
%       YouTube: https://www.youtube.com/watch?v=T6BoYRvQ6rg
%
% -------------------------------------------------------------------------
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   date: 2021-04-17
%   modified: 2022-02-16

% print status
time_start = tic;
str_date_time = sprintf( '%04d-%02d-%02d: %02d:%02d:%02d', fix( clock ) );
fprintf( '\t %s: Delay-and-Sum (DAS) Beamforming [ Fourier domain, steered plane wave ]... ', str_date_time );

%--------------------------------------------------------------------------
% 1.) check arguments
%--------------------------------------------------------------------------
% ensure at least 9 and at most 13 arguments
narginchk( 9, 13 );

% ensure existence of nonempty index_t0
if nargin < 10 || isempty( index_t0 )
	index_t0 = 0;
end

% ensure nonnegative integer for index_t0
mustBeNonnegative( index_t0 );
mustBeInteger( index_t0 );
mustBeScalarOrEmpty( index_t0 );

% ensure existence of nonempty window
if nargin < 11 || isempty( window )
    window = windows.tukey( 0.2 );
end

% ensure class windows.window
if ~isa( window, 'windows.window' )
	errorStruct.message = 'window must be windows.window!';
	errorStruct.identifier = 'das_pw:NoWindow';
	error( errorStruct );
end

% ensure existence of nonempty F_number
if nargin < 12 || isempty( F_number )
    F_number = f_numbers.grating.angle_lb( 60, 1.5 );
end

% ensure class f_numbers.f_number
if ~isa( F_number, 'f_numbers.f_number' )
	errorStruct.message = 'F_number must be f_numbers.f_number!';
	errorStruct.identifier = 'das_pw:NoFNumber';
	error( errorStruct );
end

% ensure existence of nonempty normalization
if nargin < 13 || isempty( normalization )
    normalization = normalizations.on;
end

% ensure class normalizations.normalization
if ~isa( normalization, 'normalizations.normalization' )
    errorStruct.message = 'normalization must be normalizations.normalization!';
    errorStruct.identifier = 'das_pw:NoNormalization';
    error( errorStruct );
end

%--------------------------------------------------------------------------
% 2.) geometry
%--------------------------------------------------------------------------
% number of array elements
N_elements = size( data_RF, 2 );
M_elements = ( N_elements - 1 ) / 2;

% half-width of the elements
element_width_over_two = element_width / 2;

% centroids of element faces
positions_ctr_x = (-M_elements:M_elements) * element_pitch;
positions_lbs_x = positions_ctr_x - element_width_over_two;
positions_ubs_x = positions_ctr_x + element_width_over_two;

% reference position
position_ctr_x_ref = sign( e_theta( 1 ) ) * positions_ctr_x( 1 );

% lateral distances [ N_pos_x, N_elements ]
dist_lateral = positions_ctr_x - positions_x( : );
indicator_lower = dist_lateral < 0;

% incident wave travel times
t_in_plus_shift = ( e_theta( 1 ) * ( positions_x - position_ctr_x_ref ) + e_theta( 2 ) * positions_z( : ) ) / c_0 + index_t0 / f_s;

% scattered wave travel times
t_sc_lateral_squared = ( dist_lateral / c_0 ).^2;
t_sc_axial_squared = ( positions_z / c_0 ).^2;

% maximum relative time shift in the electronic focusing ( zero padding in DFT )
% N_samples_t_add = ceil( ( sqrt( ( 2 * M_elements * element_pitch )^2 + positions_z( 1 )^2 ) - positions_z( 1 ) ) * f_s / c_0 );
% effect of bandpass filter: ceil( 10^( 60 / 20 ) / ( pi * diff( f_bounds ) ) * f_s )
N_samples_t_add = 2000;

%--------------------------------------------------------------------------
% 3.) create frequency axis
%--------------------------------------------------------------------------
% number of samples in time domain
N_samples_t = size( data_RF, 1 );

% ensure odd number of points in DFT
N_points_dft = N_samples_t + N_samples_t_add + 1 - mod( N_samples_t + N_samples_t_add, 2 );

% boundary frequency indices
index_Omega_lb = ceil( f_bounds( 1 ) * N_points_dft / f_s ) + 1;
index_Omega_ub = floor( f_bounds( 2 ) * N_points_dft / f_s ) + 1;
indices_Omega = (index_Omega_lb:index_Omega_ub).';

% frequency axes
axis_f_bp = f_s * ( indices_Omega - 1 ) / N_points_dft;
axis_omega_bp = 2 * pi * axis_f_bp;

% number of relevant discrete frequencies
N_samples_f = numel( indices_Omega );

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
% 5.) electronic focusing
%--------------------------------------------------------------------------
% compute DFT of RF data
data_RF_dft = fft( data_RF, N_points_dft, 1 );
data_RF_dft_analy = 2 * data_RF_dft( indices_Omega, : );

% desired half-widths of the receive subapertures
width_aperture_over_two_desired = positions_z ./ ( 2 * F_number_values );

% initialize image w/ zeros
image = zeros( numel( positions_z ), numel( positions_x ) );

window_f = windows.tukey( 0.2 );

% iterate lateral voxel positions
for index_pos_x = 1:numel( positions_x )

    % print progress in percent
    fprintf( '%5.1f %%', ( index_pos_x - 1 ) / numel( positions_x ) * 1e2 );

    % map desired bounds of the receive aperture to element indices
    indices_aperture_lb = max( ceil( M_elements + ( positions_x( index_pos_x ) - width_aperture_over_two_desired ) / element_pitch ) + 1, 1 );
    indices_aperture_ub = min( floor( M_elements + ( positions_x( index_pos_x ) + width_aperture_over_two_desired ) / element_pitch ) + 1, N_elements );

    % actual width of the receive aperture
    width_aperture_lower_over_two = positions_x( index_pos_x ) - positions_lbs_x( indices_aperture_lb );
    width_aperture_upper_over_two = positions_ubs_x( indices_aperture_ub ) - positions_x( index_pos_x );

    % iterate axial voxel positions
    for index_pos_z = 1:numel( positions_z )

        %------------------------------------------------------------------
        % a) compute time-of-flight (TOF)
        %------------------------------------------------------------------
        t_sc = sqrt( t_sc_lateral_squared( index_pos_x, : ) + t_sc_axial_squared( index_pos_z ) );
        tof = t_in_plus_shift( index_pos_z, index_pos_x ) + t_sc;

        weights = exp( 1j * axis_omega_bp * tof );

        %------------------------------------------------------------------
        % b) compute apodization weights for current voxel
        %------------------------------------------------------------------
        % check type of window function
        if isa( window, 'windows.boxcar' )

            %--------------------------------------------------------------
            % i.) simple solution for boxcar window
            %--------------------------------------------------------------
            window_samples = compute_samples( window, dist_lateral( index_pos_x, : ) ./ width_aperture_over_two_desired( :, index_pos_z ) );

        else

            %--------------------------------------------------------------
            % ii.) complex solution for other windows
            %--------------------------------------------------------------
            % sample window functions
            indicator_lower_act = indicator_lower( index_pos_x, : );
            window_samples_lower = compute_samples( window, dist_lateral( index_pos_x, indicator_lower_act ) ./ width_aperture_lower_over_two( :, index_pos_z ) );
            window_samples_upper = compute_samples( window, dist_lateral( index_pos_x, ~indicator_lower_act ) ./ width_aperture_upper_over_two( :, index_pos_z ) );
            window_samples = [ window_samples_lower, window_samples_upper ];

        end % if isa( window, 'windows.boxcar' )

        % normalize window samples
        window_samples = apply( normalization, window_samples );

        % apply window
%         lengths_over_two = sum( window_samples > 0, 1 ) / 2;
%         temp = compute_samples( window_f, ( (0:N_samples_f-1).' - lengths_over_two ) ./ lengths_over_two );
%         window_samples = window_samples .* temp;

        %------------------------------------------------------------------
        % c) apply weights and focus RF signals
        %------------------------------------------------------------------
        data_RF_dft_analy_focused = sum( window_samples .* weights .* data_RF_dft_analy, 2 );

%         figure( index_pos_z + 3 );
%         temp = zeros( index_Omega_ub, N_elements );
%         temp( indices_Omega, : ) = window_samples .* weights;
%         window_samples_td = ifft( temp, N_points_dft, 1, 'symmetric' );
%         imagesc( illustration.dB( window_samples_td, 20 ), [ -60, 0 ] );
%         temp = zeros( index_Omega_ub, 1 );
%         temp( indices_Omega ) = data_RF_dft_analy_focused;
%         temp = ifft( temp, N_points_dft, 1, 'symmetric' );
%         plot( (1:N_samples_t), temp, (1:N_samples_t), abs( hilbert( temp ) ) );

        %------------------------------------------------------------------
        % d) inverse DFT of focused signal at time index 0
        %------------------------------------------------------------------
        image( index_pos_z, index_pos_x ) = sum( data_RF_dft_analy_focused, 1 );

	end % for index_pos_z = 1:numel( positions_z )

	% erase progress in percent
	fprintf( '\b\b\b\b\b\b\b' );

end % for index_pos_x = 1:numel( positions_x )

% return focused RF signal
if nargout >= 3
	signal = zeros( N_points_dft, 1 );
	signal( indices_Omega ) = data_RF_dft_analy_focused;
    signal = ifft( signal, N_points_dft, 1 );
end

% infer and print elapsed time
time_elapsed = toc( time_start );
fprintf( 'done! (%f s)\n', time_elapsed );

end % function [ image, F_number_values, signal ] = das_pw( positions_x, positions_z, data_RF, f_s, e_theta, element_width, element_pitch, f_bounds, c_0, index_t0, window, F_number, normalization )
