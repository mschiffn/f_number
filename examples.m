% exemplary applications of
% the frequency-dependent F-number to
% coherent plane-wave compounding
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
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   date: 2023-07-14
%   modified: 2023-08-02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0.) parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load RF data (tissue phantom)
load( 'data_RF.mat' );

% specify bandwidth
f_bounds = [ f_lb, f_ub ];

% specify time index of the sample extracted from the focused RF signal
index_t0 = 8;

% steering angles in rad
steering_angles_rad = deg2rad( steering_angles_deg );

% dependent parameters
positions_x = (-300:300) * element_pitch / 4;
positions_z = ( 64 + (0:511) ) * element_pitch / 4;

% dynamic range for all illustrations
dynamic_range_dB = 60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.) typical usage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% typical usage of
% the Fourier-domain beamformer in
% coherent plane-wave compounding

%--------------------------------------------------------------------------
% 1.) compute images
%--------------------------------------------------------------------------
% specify cell array for results
images = cell( 1, numel( steering_angles_rad ) );

% iterate steering angles
for index_angle = 1:numel( steering_angles_rad )

    % call Fourier-domain beamformer
    images{ index_angle } = das_pw( positions_x, positions_z, data_RF( :, :, index_angle ), f_s, steering_angles_rad( index_angle ), element_width, element_pitch, c_avg, f_bounds, index_t0 );

end

% compute compound image
image_compound = sum( cat( 3, images{ : } ), 3 );

%--------------------------------------------------------------------------
% 2.) show results
%--------------------------------------------------------------------------
c_limits = [ -dynamic_range_dB, 0 ];

figure( 1 );
for index_angle = 1:numel( steering_angles_rad )
    subplot( 1, numel( steering_angles_rad ), index_angle );
    imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( images{ index_angle } ) / max( abs( images{ index_angle }( : ) ) ) ), c_limits );
    title( sprintf( 'Steering angle: %.0f°', steering_angles_deg( index_angle ) ) );
    xlabel( 'Lateral position (mm)' );
    ylabel( 'Axial position (mm)' );
    colormap gray;
end

figure( 2 );
imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( image_compound ) / max( abs( image_compound( : ) ) ) ), c_limits );
title( 'Coherent compounding' );
xlabel( 'Lateral position (mm)' );
ylabel( 'Axial position (mm)' );
colormap gray;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.) graphical abstract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproduce the graphical abstract of [1].

% select steering angle of -20°
index_angle = 1;

% specify receive F-numbers
F_numbers_rx = { f_numbers.constant( 0 ), f_numbers.constant( 1.5 ), f_numbers.grating.angle_lb( 60, 1.5 ) };
str_titles = { 'Full aperture', 'Fixed', 'Proposed' };

%--------------------------------------------------------------------------
% 1.) compute images
%--------------------------------------------------------------------------
% specify cell array for results
images = cell( 1, numel( F_numbers_rx ) );

% iterate F-numbers
for index_F = 1:numel( F_numbers_rx )

    % call Fourier-domain beamformer
    images{ index_F } = das_pw( positions_x, positions_z, data_RF( :, :, index_angle ), f_s, steering_angles_rad( index_angle ), element_width, element_pitch, c_avg, f_bounds, index_t0, [], F_numbers_rx{ index_F } );

end

%--------------------------------------------------------------------------
% 2.) show results
%--------------------------------------------------------------------------
c_limits = [ -dynamic_range_dB, 0 ];

figure( 2 );
for index_F = 1:numel( F_numbers_rx )
    subplot( 1, numel( F_numbers_rx ), index_F );
    imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( images{ index_F } ) / max( abs( images{ index_F }( : ) ) ) ), c_limits );
    title( str_titles{ index_F }, 'Interpreter', 'none' );
    xlabel( 'Lateral position (mm)' );
    ylabel( 'Axial position (mm)' );
    colormap gray;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.) apodization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select steering angle of -20°
index_angle = 1;

% specify window functions
windows_rx = { windows.boxcar, windows.tukey( 0.2 ), windows.triangular, windows.hann };

%--------------------------------------------------------------------------
% 1.) compute images
%--------------------------------------------------------------------------
% specify cell array for results
images = cell( 1, numel( windows_rx ) );

% iterate window functions
for index_window = 1:numel( windows_rx )

    % call Fourier-domain beamformer
    images{ index_window } = das_pw( positions_x, positions_z, data_RF( :, :, index_angle ), f_s, steering_angles_rad( index_angle ), element_width, element_pitch, c_avg, f_bounds, index_t0, windows_rx{ index_window } );

end

%--------------------------------------------------------------------------
% 2.) show results
%--------------------------------------------------------------------------
c_limits = [ -dynamic_range_dB, 0 ];

figure( 3 );
for index_window = 1:numel( windows_rx )
    subplot( 1, numel( windows_rx ), index_window );
    imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( images{ index_window } ) / max( abs( images{ index_window }( : ) ) ) ), c_limits );
    title( sprintf( '%s', windows_rx{ index_window } ), 'Interpreter', 'none' );
    xlabel( 'Lateral position (mm)' );
    ylabel( 'Axial position (mm)' );
    colormap gray;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4.) normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select steering angle of -20°
index_angle = 1;

% specify types of normalization
normalization_rx = { normalizations.off, normalizations.on };

%--------------------------------------------------------------------------
% 1.) compute images
%--------------------------------------------------------------------------
% specify cell array for results
images = cell( 1, numel( normalization_rx ) );

% iterate window functions
for index_normalization = 1:numel( normalization_rx )

    % call Fourier-domain beamformer
    images{ index_normalization } = das_pw( positions_x, positions_z, data_RF( :, :, index_angle ), f_s, steering_angles_rad( index_angle ), element_width, element_pitch, c_avg, f_bounds, index_t0, [], [], normalization_rx{ index_normalization } );

end

%--------------------------------------------------------------------------
% 2.) show results
%--------------------------------------------------------------------------
c_limits = [ -dynamic_range_dB, 0 ];

figure( 4 );
for index_normalization = 1:numel( normalization_rx )
    subplot( 1, numel( normalization_rx ), index_normalization );
    imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( images{ index_normalization } ) / max( abs( images{ index_normalization }( : ) ) ) ), c_limits );
    title( sprintf( 'Normalization: %s', normalization_rx{ index_normalization } ), 'Interpreter', 'none' );
    xlabel( 'Lateral position (mm)' );
    ylabel( 'Axial position (mm)' );
    colormap gray;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5.) bandwidth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The frequency bounds specify
% the bandwidth for
% all computations.
%
% Smaller bandwidths reduce
% the computational costs at
% the expense of
% the image quality.

% select steering angle of -20°
index_angle = 1;

% bandwidths used for image formation
f_c = mean( f_bounds );
fractional_bandwidths = [ 1, 0.75, 0.5, 0.25 ];

%--------------------------------------------------------------------------
% 1.) compute images
%--------------------------------------------------------------------------
% specify cell array for results
images = cell( 1, numel( fractional_bandwidths ) );

% iterate fractional bandwidths
for index_bandwidth = 1:numel( fractional_bandwidths )

    % current frequency bounds
    f_bounds_act = f_c * ( 1 + [ -1, 1 ] * fractional_bandwidths( index_bandwidth ) / 2 );

    % call Fourier-domain beamformer
    images{ index_bandwidth } = das_pw( positions_x, positions_z, data_RF( :, :, index_angle ), f_s, steering_angles_rad( index_angle ), element_width, element_pitch, c_avg, f_bounds_act, index_t0 );

end

%--------------------------------------------------------------------------
% 2.) show results
%--------------------------------------------------------------------------
c_limits = [ -dynamic_range_dB, 0 ];

figure( 5 );
for index_bandwidth = 1:numel( fractional_bandwidths )
    subplot( 2, numel( fractional_bandwidths ), index_bandwidth );
    imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( images{ index_bandwidth } ) / max( abs( images{ index_bandwidth }( : ) ) ) ), c_limits );
    title( sprintf( 'Fractional bandwidth: %d %%', fractional_bandwidths( index_bandwidth ) * 1e2 ) );
    xlabel( 'Lateral position (mm)' );
    ylabel( 'Axial position (mm)' );
    colormap gray;
    subplot( 2, numel( fractional_bandwidths ), index_bandwidth + numel( fractional_bandwidths ) );
    imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( images{ index_bandwidth } ) / max( abs( images{ index_bandwidth }( : ) ) ) ), c_limits );
end
