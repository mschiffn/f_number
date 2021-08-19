function [ image, F_number_values, signal ] = das_pw( positions_x, positions_z, data_RF, f_s, e_theta, element_width, element_pitch, M_elements_start, f_bounds, c_0, index_t0, window, F_number )
%DAS_PW Delay-and-Sum (DAS) Beamforming [ Fourier domain, steered plane wave ]
%
% Computes a sonogram using the delay-and-sum (DAS) algorithm in
% the Fourier domain.
% This domain enables the exact application of
% all time shifts independent of
% the sampling rate.
% A uniform linear transducer array is assumed.
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% 01.) positions_x:         lateral pixel positions (m)
% 02.) positions_z:         axial pixel positions (m)
% 03.) data_RF:             RF data (2d array; 1st dimension: time, 2nd dimension: array element index)
% 04.) f_s:                 sampling rate of the RF data (Hz)
% 05.) e_theta:             propagation direction of incident plane wave (1)
% 06.) element_width:       element width of the linear transducer array (m)
% 07.) element_pitch:       element pitch of the linear transducer array (m)
% 08.) M_elements_start:	index of the first element (1)
% 09.) f_bounds:            frequency bounds (Hz)
% 10.) c_0:                 speed of sound (m/s)
% 11.) index_t0:            time index of the sample extracted from the focused RF signal (1)
% 12.) window:              window function for receive apodization (  object of class windows.window )
% 13.) F_number:            receive F-number (  object of class f_numbers.f_number )
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% 01.) image:               complex-valued DAS image
% 02.) F_number_values:     computed F-number for each frequency
% 03.) signal:              focused RF signal for last image voxel
%
% REFERENCES:
%	[1] M. F. Schiffner and G. Schmitz, “Frequency-Dependent F-Number Increases the Contrast and the Spatial Resolution,” in
%       2021 IEEE Int. Ultrasonics Symp. (IUS), 2021.
%       DOI:
%
% ABOUT:
%	author: Martin F. Schiffner
%	date: 2021-04-17
%	modified: 2021-08-07

% TODO: special case: constant F-number
% TODO: maximum half-width of the aperture?

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
if nargin < 11 || isempty( index_t0 )
	index_t0 = 0;
end

% ensure nonnegative integer for index_t0
mustBeNonnegative( index_t0 );
mustBeInteger( index_t0 );
mustBeScalarOrEmpty( index_t0 );

% ensure existence of nonempty window
if nargin < 12 || isempty( window )
	window = windows.boxcar;
end

% ensure class windows.window
if ~isa( window, 'windows.window' )
	errorStruct.message = 'window must be windows.window!';
	errorStruct.identifier = 'das_pw:NoWindow';
	error( errorStruct );
end

% ensure existence of nonempty F_number
if nargin < 13 || isempty( F_number )
	F_number = f_numbers.grating.angle_lb( 60, 1 );
end

% ensure class f_numbers.f_number
if ~isa( F_number, 'f_numbers.f_number' )
	errorStruct.message = 'F_number must be f_numbers.f_number!';
	errorStruct.identifier = 'das_pw:NoFNumber';
	error( errorStruct );
end

%--------------------------------------------------------------------------
% 2.) geometry
%--------------------------------------------------------------------------
% centroids of element faces
positions_ctr_x = (M_elements_start:-M_elements_start) * element_pitch;

% reference position
if e_theta( 1 ) >= 0
	position_ctr_x_ref = positions_ctr_x( 1 );
else
	position_ctr_x_ref = positions_ctr_x( end );
end

% lateral distances
dist_lateral = positions_ctr_x - positions_x( : );

% incident wave travel times
t_in = ( e_theta( 1 ) * ( positions_x - position_ctr_x_ref ) + e_theta( 2 ) * positions_z( : ) ) / c_0;

% scattered wave travel times
t_sc_lateral_squared = ( dist_lateral / c_0 ).^2;
t_sc_axial_squared = ( positions_z / c_0 ).^2;

% maximum relative time shift in the electronic focusing ( zero padding in DFT )
% N_samples_t_add = ceil( ( sqrt( ( 2 * abs( M_elements_start ) * element_pitch )^2 + positions_z( 1 )^2 ) - positions_z( 1 ) ) * f_s / c_0 );
N_samples_t_add = 0;

% permissible maximum half-width of the receive subaperture
width_aperture_over_2_max = 2 * abs( M_elements_start ) * element_pitch + element_width;
% width_aperture_over_2_min = ( 3 * element_pitch + element_width ) / 2;
% width_aperture_over_2_max = abs( M_elements_start ) * element_pitch + 0.5 * element_width;

%--------------------------------------------------------------------------
% 3.) create frequency axis
%--------------------------------------------------------------------------
% number of samples in time domain
N_samples_t = size( data_RF, 1 );

% ensure odd number of points in DFT
N_points_dft = N_samples_t + N_samples_t_add + 1 - mod( N_samples_t + N_samples_t_add, 2 );

% boundary frequency indices
index_Omega_lb = ceil( f_bounds( 1 ) * N_points_dft / f_s );
index_Omega_ub = floor( f_bounds( 2 ) * N_points_dft / f_s );
indices_Omega = (index_Omega_lb:index_Omega_ub).';

% frequency axes
axis_f_bp = f_s * indices_Omega / N_points_dft;
axis_omega_bp = 2 * pi * axis_f_bp;

%--------------------------------------------------------------------------
% 4.) frequency-dependent F-number
%--------------------------------------------------------------------------
% element pitch-to-wavelength ratio
element_pitch_over_lambda = element_pitch * axis_f_bp / c_0;

% compute values of the F-number for each frequency
F_number_values = compute_values( F_number, element_pitch_over_lambda );

%--------------------------------------------------------------------------
% 5.) electronic focusing
%--------------------------------------------------------------------------
% compute DFT of RF data
data_RF_dft = fft( data_RF, N_points_dft, 1 );
data_RF_dft_analy = 2 * data_RF_dft( index_Omega_lb:index_Omega_ub, : );

% half-widths of the receive subapertures
width_aperture_over_2 = min( width_aperture_over_2_max, positions_z ./ ( 2 * F_number_values ) );

% initialize image w/ zeros
image = zeros( numel( positions_z ), numel( positions_x ) );

% iterate lateral voxel positions
for index_pos_x = 1:numel( positions_x )

	% print progress in percent
	fprintf( '%5.1f %%', ( index_pos_x - 1 ) / numel( positions_x ) * 1e2 );

	% iterate axial voxel positions
	for index_pos_z = 1:numel( positions_z )

        %------------------------------------------------------------------
        % a) compute time-of-flight (TOF)
        %------------------------------------------------------------------
        t_sc = sqrt( t_sc_lateral_squared( index_pos_x, : ) + t_sc_axial_squared( index_pos_z ) );
        tof = t_in( index_pos_z, index_pos_x ) + t_sc;

        weights = exp( 1j * axis_omega_bp * ( tof + index_t0 / f_s ) );

        %------------------------------------------------------------------
        % b) compute apodization weights for current pixel
        %------------------------------------------------------------------
        % frequency-dependent apodization weights
        if isa( F_number, 'f_numbers.constant' )
            %--------------------------------------------------------------
            % constant F-number
            %--------------------------------------------------------------
            apodization_weights = compute_samples( window, dist_lateral( index_pos_x, : ), width_aperture_over_2( 1, index_pos_z ) );
        else
            %--------------------------------------------------------------
            % frequency-dependent F-number
            %--------------------------------------------------------------
            % 1. dimension: frequency
            % 2. dimension: element index
            apodization_weights = compute_samples( window, dist_lateral( index_pos_x, : ), width_aperture_over_2( :, index_pos_z ) );
        end

        apodization_weights_sum = sum( apodization_weights, 2 );
        indicator = apodization_weights_sum < eps;
        apodization_weights = apodization_weights ./ apodization_weights_sum;
        apodization_weights( indicator, : ) = 0;

        %------------------------------------------------------------------
        % c) apply weights and focus RF signals
        %------------------------------------------------------------------
        data_RF_dft_analy_focused = sum( apodization_weights .* weights .* data_RF_dft_analy, 2 );

%         figure( index_pos_z + 1 );
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

figure(998);
imagesc( illustration.dB( image, 20 ), [ -70, 0 ] );
title( string( F_number ), 'Interpreter', 'none' );
colormap gray;

% return focused RF signal
if nargout >= 3
	signal = zeros( N_points_dft, 1 );
	signal( indices_Omega ) = data_RF_dft_analy_focused;
    signal = ifft( signal, N_points_dft, 1 );
end

% infer and print elapsed time
time_elapsed = toc( time_start );
fprintf( 'done! (%f s)\n', time_elapsed );

end % function [ image, F_number_values, signal ] = das_pw( positions_x, positions_z, data_RF, f_s, e_theta, element_width, element_pitch, M_elements_start, f_bounds, c_0, index_t0, apodization, F_number )
