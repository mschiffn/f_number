function [ images, image_compound, image_compound_linear, image_compound_linear_normalized ] = das_prog_scan( pos_x_focus_tx, pos_z_focus_rx, data_RF, f_s, pos_z_focus_tx, element_width, element_pitch, c_ref, f_bounds, index_t0, window_tx, window_rx, F_number_tx, F_number_rx, normalization_tx, normalization_rx, dyn_rx_focus )
% das_prog_scan Synthesize the RF data acquired by multifocal progressive scanning
%
% This function synthesizes
% the RF data recorded during
% multiple monofocal progressive scans from
% a complete synthetic aperture scan.
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% REQUIRED
%   01.) pos_x_focus_tx:    lateral positions of the transmit foci (m)
%   02.) pos_z_focus_rx:    axial positions of the receive foci (m)
%   03.) data_RF:           RF data of a full SA scan (3d array; 1st dimension: time, 2nd dimension: rx array element index, 3rd dimension: tx array element index)
%   04.) f_s:               sampling rate of the RF data (Hz)
%   05.) pos_z_focus_tx:    axial positions of the transmit foci (m)
%   06.) element_width:     element width of the uniform linear array (m)
%   07.) element_pitch:     element pitch of the uniform linear array (m)
%   08.) c_ref:             speed of sound (m/s)
%
% OPTIONAL
%   09.) f_bounds:          frequency bounds (Hz)
%   10.) index_t0:          time index of the sample extracted from the focused RF signal (1)
%   11.) window_tx:         window function for transmit apodization ( object of class windows.window )
%   12.) window_rx:         window function for receive apodization ( object of class windows.window )
%   13.) F_number_tx:       transmit F-number ( object of class f_numbers.f_number )
%   14.) F_number_rx:       receive F-number ( object of class f_numbers.f_number )
%   15.) normalization_tx:  normalization method for transmit apodization weights ( object of class normalizations.normalization )
%   16.) normalization_rx:  normalization method for receive apodization weights ( object of class normalizations.normalization )
%   17.) dyn_rx_focus:      toggle dynamic receive focusing (default: true)
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%   01.) images:          complex-valued and focused RF lines (cell array; one entry for each transmit focal length; 1st dimension: time, 2nd dimension: lateral position)
%   02.) image_compound:                    compounded DAS image [ no weighting ]
%   03.) image_compound_linear:             compounded DAS image [ linear weighting ]
%   04.) image_compound_linear_normalized:  compounded DAS image [ linear weighting, normalized ]
%
% -------------------------------------------------------------------------
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   created: 2012-09-13
%   modified: 2025-04-03

%--------------------------------------------------------------------------
% 1.) check arguments
%--------------------------------------------------------------------------
% ensure valid number of arguments
narginchk( 8, 17 );

N_WORKERS = 6;

% ensure existence of nonempty f_bounds
if nargin < 9 || isempty( f_bounds )
    f_bounds = [ eps, f_s / 2 ];
end

% ensure existence of nonempty index_t0
if nargin < 10 || isempty( index_t0 )
    index_t0 = 0;
end

% ensure existence of nonempty tx window
if nargin < 11 || isempty( window_tx )
    window_tx = windows.tukey( 0.2 );
end

% ensure class windows.window
if ~( isa( window_tx, 'windows.window' ) && isscalar( window_tx ) )
    errorStruct.message = 'window_tx must be scalar windows.window!';
    errorStruct.identifier = 'synthesize_progressive_scan:NoTXWindow';
    error( errorStruct );
end

% ensure existence of nonempty rx window
if nargin < 12 || isempty( window_rx )
    window_rx = windows.tukey( 0.2 );
end

% ensure class windows.window
if ~( isa( window_rx, 'windows.window' ) && isscalar( window_rx ) )
    errorStruct.message = 'window_rx must be scalar windows.window!';
    errorStruct.identifier = 'synthesize_progressive_scan:NoRXWindow';
    error( errorStruct );
end

% ensure existence of nonempty tx F_number
if nargin < 13 || isempty( F_number_tx )
    F_number_tx = f_numbers.grating.angle_lb( 60, 1.5 );
end

% ensure class f_numbers.f_number
if ~( isa( F_number_tx, 'f_numbers.f_number' ) && isscalar( F_number_tx ) )
	errorStruct.message = 'F_number_tx must be scalar f_numbers.f_number!';
	errorStruct.identifier = 'synthesize_progressive_scan:NoTXFNumber';
	error( errorStruct );
end

% ensure existence of nonempty rx F_number
if nargin < 14 || isempty( F_number_rx )
    F_number_rx = f_numbers.grating.angle_lb( 60, 1.5 );
end

% ensure class f_numbers.f_number
if ~( isa( F_number_rx, 'f_numbers.f_number' ) && isscalar( F_number_rx ) )
	errorStruct.message = 'F_number_rx must be scalar f_numbers.f_number!';
	errorStruct.identifier = 'synthesize_progressive_scan:NoRXFNumber';
	error( errorStruct );
end

% ensure existence of nonempty tx normalization
if nargin < 15 || isempty( normalization_tx )
    normalization_tx = normalizations.on;
end

% ensure class normalizations.normalization
if ~( isa( normalization_tx, 'normalizations.normalization' ) && isscalar( normalization_tx ) )
    errorStruct.message = 'normalization_tx must be scalar normalizations.normalization!';
    errorStruct.identifier = 'synthesize_progressive_scan:NoTXNormalization';
    error( errorStruct );
end

% ensure existence of nonempty rx normalization
if nargin < 16 || isempty( normalization_rx )
    normalization_rx = normalizations.on;
end

% ensure class normalizations.normalization
if ~( isa( normalization_rx, 'normalizations.normalization' ) && isscalar( normalization_rx ) )
    errorStruct.message = 'normalization_rx must be scalar normalizations.normalization!';
    errorStruct.identifier = 'synthesize_progressive_scan:NoRXNormalization';
    error( errorStruct );
end

% ensure nonempty dyn_rx_focus
if nargin < 17 || isempty( dyn_rx_focus )
    dyn_rx_focus = true;
end

%--------------------------------------------------------------------------
% 2.) geometrical parameters
%--------------------------------------------------------------------------
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
dist_lateral = positions_ctr_x - pos_x_focus_tx.';

% determine left part of the aperture
indicator_left = dist_lateral < 0;

% number of lateral and axial focus positions (tx)
N_focus_x_tx = numel( pos_x_focus_tx );     % number of lines per image
N_focus_z_tx = numel( pos_z_focus_tx );     % number of images

% number of axial focus positions (rx)
N_pos_z_rx = numel( pos_z_focus_rx );       % number of samples per line

%--------------------------------------------------------------------------
% 3.) order of DFT, temporal frequency axes, bandpass filter
%--------------------------------------------------------------------------
N_points_dft = size( data_RF, 1 );
% employ zero-padding to get odd signal length
if mod( N_points_dft, 2 ) == 0
    N_points_dft = N_points_dft + 1;
end
% assertion: N_points_dft is odd

% boundary frequency indices
index_f_lb = ceil( f_bounds( 1 ) * N_points_dft / f_s ) + 1;
index_f_ub = floor( f_bounds( 2 ) * N_points_dft / f_s ) + 1;
indices_f = (index_f_lb:index_f_ub).';

% frequency axes
axis_f_bp = f_s * ( indices_f - 1 ) / N_points_dft;
axis_Omega_bp = 2 * pi * ( indices_f - 1 ) / N_points_dft;

% compute DFTs (employ zero-padding)
data_RF_SA_dft = fft( data_RF, N_points_dft, 1 );

% compute DFTs of analytic signals in the relevant frequency band
data_RF_SA_analy_dft = 2 * data_RF_SA_dft( indices_f, :, : );

%--------------------------------------------------------------------------
% 4.) dynamic aperture
%--------------------------------------------------------------------------
% element pitch-to-wavelength ratio
element_pitch_over_lambda = element_pitch * axis_f_bp / c_ref;

% compute values of the F-numbers for each frequency
F_number_tx_values = compute_values( F_number_tx, element_pitch_over_lambda );
F_number_rx_values = compute_values( F_number_rx, element_pitch_over_lambda );

% maximum F-numbers for each axial position
F_tx_ub = sqrt( axis_f_bp .* pos_z_focus_tx / ( 2.88 * c_ref ) );
F_rx_ub = sqrt( axis_f_bp .* pos_z_focus_rx / ( 2.88 * c_ref ) );

% desired half-widths of the apertures
width_aperture_tx_over_two_desired = pos_z_focus_tx ./ ( 2 * F_number_tx_values );
width_aperture_rx_over_two_desired = pos_z_focus_rx ./ ( 2 * F_number_rx_values );

%--------------------------------------------------------------------------
% 5.) compute progressive scan
%--------------------------------------------------------------------------
% specify cell arrays for images
images = cell( 1, N_focus_z_tx );
images_normalized = cell( 1, N_focus_z_tx );

% parpool( N_WORKERS );

% iterate axial tx foci
% parfor index_focus_tx_z = 1:N_focus_z_tx
for index_focus_tx_z = 1:N_focus_z_tx

    fprintf( 'index_focus_tx_z = %d of %d\n', index_focus_tx_z, N_focus_z_tx );

    % allocate memory for current image
    data_RF_act = zeros( N_pos_z_rx, N_focus_x_tx );

    % iterate lateral tx foci
    for index_focus_tx_x = 1:N_focus_x_tx

        fprintf( '\t index_focus_tx_x = %d of %d\n', index_focus_tx_x, N_focus_x_tx );

        %------------------------------------------------------------------
        % a) transmit focusing (use all elements for rx)
        %------------------------------------------------------------------
        % map desired bounds of the tx aperture to element indices
        indices_aperture_tx_lb = max( ceil( M_elements + ( pos_x_focus_tx( index_focus_tx_x ) - width_aperture_tx_over_two_desired( :, index_focus_tx_z ) ) / element_pitch ) + 1, 1 );
        indices_aperture_tx_ub = min( floor( M_elements + ( pos_x_focus_tx( index_focus_tx_x ) + width_aperture_tx_over_two_desired( :, index_focus_tx_z ) ) / element_pitch ) + 1, N_elements );

        % frequency-dependent processing of each z-coordinate
        indicator_tx_valid = indices_aperture_tx_lb <= indices_aperture_tx_ub;

        % skip iteration for inactive tx aperture
        if ~any( indicator_tx_valid ), continue, end

        % actual width of the tx aperture
        width_aperture_tx_left_over_two = pos_x_focus_tx( index_focus_tx_x ) - positions_lbs_x( indices_aperture_tx_lb( indicator_tx_valid ) );
        width_aperture_tx_right_over_two = positions_ubs_x( indices_aperture_tx_ub( indicator_tx_valid ) ) - pos_x_focus_tx( index_focus_tx_x );

        % sample window functions
        indicator_left_act = indicator_left( index_focus_tx_x, : );
        window_samples_left = compute_samples( window_tx, dist_lateral( index_focus_tx_x, indicator_left_act ) ./ width_aperture_tx_left_over_two );
        window_samples_right = compute_samples( window_tx, dist_lateral( index_focus_tx_x, ~indicator_left_act ) ./ width_aperture_tx_right_over_two );
        window_samples = [ window_samples_left, window_samples_right ];

        % normalize window samples
        window_samples = reshape( apply( normalization_tx, window_samples ), [ sum( indicator_tx_valid ), 1, N_elements ] );

        % time shifts for active tx elements
        distances_tx = sqrt( ( positions_ctr_x - pos_x_focus_tx( index_focus_tx_x ) ).^2 + pos_z_focus_tx( index_focus_tx_z )^2 );
        distances_tx_min = min( distances_tx );
        N_samples_shift_focus = f_s * ( distances_tx_min - distances_tx ) / c_ref;

% TODO: why acausal?
        % apply delays for tx (acausal) and summation
        phase_shift = exp( -1j * axis_Omega_bp( indicator_tx_valid ) .* reshape( N_samples_shift_focus, [ 1, 1, N_elements ] ) );
        data_RF_analy_tx_dft = sum( window_samples .* phase_shift .* data_RF_SA_analy_dft( indicator_tx_valid, :, : ), 3 );

        %------------------------------------------------------------------
        % b) receive focusing
        %------------------------------------------------------------------
        if dyn_rx_focus == 1

            %--------------------------------------------------------------
            % dynamic receive focusing
            %--------------------------------------------------------------
            distances_rx = sqrt( repmat( ( positions_ctr_x - pos_x_focus_tx( index_focus_tx_x ) ).^2, [ N_pos_z_rx, 1 ] ) + repmat( ( pos_z_focus_rx.^2 )', [ 1, N_elements ] ) );
            distances_rx_min = min( distances_rx, [], 2 );
            N_samples_shift_focus = f_s * ( distances_rx_min - distances_rx ) / c_ref;

            % iterate axial rx foci
            for index_pos_z_rx = 1:N_pos_z_rx

                %----------------------------------------------------------
                % determine active rx aperture
                %----------------------------------------------------------
                % map desired bounds of the rx aperture to element indices
                indices_aperture_rx_lb = max( ceil( M_elements + ( pos_x_focus_tx( index_focus_tx_x ) - width_aperture_rx_over_two_desired( indicator_tx_valid, index_pos_z_rx ) ) / element_pitch ) + 1, 1 );
                indices_aperture_rx_ub = min( floor( M_elements + ( pos_x_focus_tx( index_focus_tx_x ) + width_aperture_rx_over_two_desired( indicator_tx_valid, index_pos_z_rx ) ) / element_pitch ) + 1, N_elements );

                % frequency-dependent processing of each z-coordinate
                indicator_rx_valid = indices_aperture_rx_lb <= indices_aperture_rx_ub;

                % skip iteration for inactive rx aperture
                if ~any( indicator_rx_valid ), continue, end

                % actual width of the tx aperture
                width_aperture_rx_left_over_two = pos_x_focus_tx( index_focus_tx_x ) - positions_lbs_x( indices_aperture_rx_lb( indicator_rx_valid ) );
                width_aperture_rx_right_over_two = positions_ubs_x( indices_aperture_rx_ub( indicator_rx_valid ) ) - pos_x_focus_tx( index_focus_tx_x );

                % sample window functions
                window_samples_left = compute_samples( window_rx, dist_lateral( index_focus_tx_x, indicator_left_act ) ./ width_aperture_rx_left_over_two );
                window_samples_right = compute_samples( window_rx, dist_lateral( index_focus_tx_x, ~indicator_left_act ) ./ width_aperture_rx_right_over_two );
                window_samples = [ window_samples_left, window_samples_right ];

                % normalize window samples
                window_samples = apply( normalization_rx, window_samples );

                % apply delays for rx (acausal) and summation
                phase_shift = exp( -1j * axis_Omega_bp( indicator_rx_valid ) * N_samples_shift_focus( index_pos_z_rx, : ) );
                data_RF_analy_tx_rx_dft = sum( window_samples .* phase_shift .* data_RF_analy_tx_dft( indicator_rx_valid, : ) , 2 );

                %------------------------------------------------------
                % compute single sample of analytic signal by IDFT (complex-valued!)
                %-----------------------------------------------------
                % choose correct sample to account for time shifts due to beamforming!
                index_dft = 2 * distances_rx_min( index_pos_z_rx ) * f_s / c_ref + index_t0;
                phase_shift = exp( 1j * axis_Omega_bp( indicator_rx_valid ) * ( index_dft - 1 ) );

                data_RF_act( index_pos_z_rx, index_focus_tx_x ) = sum( data_RF_analy_tx_rx_dft .* phase_shift ) / N_points_dft;

            end % for index_pos_z_rx = 1:N_pos_z_rx

        else

            %----------------------------------------------------------
            % static receive focusing
            %----------------------------------------------------------

%             distances_rx = sqrt((positions_ctr_x(index_element_rx_start:index_element_rx_stop) - pos_x_focus_tx(index_focus_tx_x)).^2 + pos_z_focus(index_focus_tx_z)^2);
%             N_samples_shift_focus = f_s * (min(distances_rx) - distances_rx) / c_ref - N_samples_shift_pulse;
%             data_RF_analy_tx_noisy_dft = data_RF_analy_tx_noisy_dft .* exp(-1j * axis_Omega_bp' * N_samples_shift_focus);
%                             
%             data_RF_act(:, index_focus_tx_x) = ifft(sum(data_RF_analy_tx_noisy_dft, 2), N_points_dft, 1, 'symmetric');
        end % if dyn_rx_focus == 1

    end % for index_focus_tx_x = 1:N_focus_x_tx

	images{ index_focus_tx_z } = data_RF_act;
	images_normalized{ index_focus_tx_z } = data_RF_act / max( abs( data_RF_act( : ) ) );

end % parfor index_focus_tx_z = 1:N_focus_z_tx

% stop worker pool
% delete(gcp('nocreate'));

%--------------------------------------------------------------------------
% 5.) create compound image (linear weighting)
%--------------------------------------------------------------------------
image_compound        = zeros( N_pos_z_rx, numel( pos_x_focus_tx ) );
image_compound_linear = zeros( N_pos_z_rx, numel( pos_x_focus_tx ) );
image_compound_linear_normalized = zeros( N_pos_z_rx, numel( pos_x_focus_tx ) );

figure(1);
hold on;
% iterate rx foci
for index_pos_z_rx = 1:N_pos_z_rx

    % find nearest tx focus
    pos_z_diff_focus_tx = pos_z_focus_rx( index_pos_z_rx ) - pos_z_focus_tx;
    [ pos_z_diff_rx_min, pos_z_diff_rx_min_index ] = min( abs( pos_z_diff_focus_tx ) );

    % compute weighting
    if pos_z_diff_focus_tx( pos_z_diff_rx_min_index ) <= 0

        %------------------------------------------------------------------
        % axial distance smaller than nearest focal point
        %------------------------------------------------------------------
        if pos_z_diff_rx_min_index > 1

            slope = 1 / ( pos_z_focus_tx( pos_z_diff_rx_min_index ) - pos_z_focus_tx( pos_z_diff_rx_min_index - 1 ) );

            weight_1 = 1 - slope * pos_z_diff_focus_tx( pos_z_diff_rx_min_index - 1 );
            weight_2 = slope * pos_z_diff_focus_tx( pos_z_diff_rx_min_index - 1 );

            fprintf( 'weight_1 = %.2f, weight_2 = %.2f, sum = %.2f, ref_1 = %d, ref_2 = %d\n', weight_1, weight_2, weight_1 + weight_2, pos_z_diff_rx_min_index - 1, pos_z_diff_rx_min_index);
            plot( index_pos_z_rx, weight_1, index_pos_z_rx, weight_2 );

            image_compound_linear( index_pos_z_rx, : ) = images{ pos_z_diff_rx_min_index - 1 }( index_pos_z_rx, : ) * weight_1 + images{ pos_z_diff_rx_min_index }( index_pos_z_rx, : ) * weight_2;
            image_compound_linear_normalized( index_pos_z_rx, : ) = images_normalized{ pos_z_diff_rx_min_index - 1 }( index_pos_z_rx, : ) * weight_1 + images_normalized{ pos_z_diff_rx_min_index }( index_pos_z_rx, : ) * weight_2;

        else

            fprintf( 'only %d\n', pos_z_diff_rx_min_index );

            plot( index_pos_z_rx, 0, index_pos_z_rx, 1 );
            % assertion: pos_z_diff_rx_min_index == 1

            image_compound_linear( index_pos_z_rx, : ) = images{ pos_z_diff_rx_min_index }( index_pos_z_rx, : );
            image_compound_linear_normalized( index_pos_z_rx, : ) = images_normalized{ pos_z_diff_rx_min_index }( index_pos_z_rx, : );

        end

    else % pos_z_diff_focus_tx(pos_z_diff_rx_min_index) > 0

        %------------------------------------------------------------------
        % axial distance larger than nearest focal point
        %------------------------------------------------------------------
        if pos_z_diff_rx_min_index < N_focus_z_tx

            slope = 1 / ( pos_z_focus_tx( pos_z_diff_rx_min_index + 1 ) - pos_z_focus_tx( pos_z_diff_rx_min_index ) );

            weight_1 = 1 - slope * pos_z_diff_focus_tx( pos_z_diff_rx_min_index );
            weight_2 = slope * pos_z_diff_focus_tx( pos_z_diff_rx_min_index );

            fprintf( 'weight_1 = %.2f, weight_2 = %.2f, sum = %.2f, ref_1 = %d, ref_2 = %d\n', weight_1, weight_2, weight_1 + weight_2, pos_z_diff_rx_min_index, pos_z_diff_rx_min_index + 1 );
            plot( index_pos_z_rx, weight_1, index_pos_z_rx, weight_2 );

            image_compound_linear( index_pos_z_rx, : ) = images{ pos_z_diff_rx_min_index }( index_pos_z_rx, : ) * weight_1 + images{ pos_z_diff_rx_min_index + 1 }( index_pos_z_rx, : ) * weight_2;
            image_compound_linear_normalized( index_pos_z_rx, : ) = images_normalized{ pos_z_diff_rx_min_index }( index_pos_z_rx, : ) * weight_1 + images_normalized{ pos_z_diff_rx_min_index + 1 }( index_pos_z_rx, : ) * weight_2;

        else

            fprintf('only %d\n', pos_z_diff_rx_min_index);
            plot(index_pos_z_rx, 1, index_pos_z_rx, 0);

            % assertion: pos_z_diff_rx_min_index == N_focus_z_tx
            image_compound_linear( index_pos_z_rx, : ) = images{ pos_z_diff_rx_min_index }( index_pos_z_rx, : );
            image_compound_linear_normalized( index_pos_z_rx, : ) = images_normalized{ pos_z_diff_rx_min_index }( index_pos_z_rx, : );

        end

    end % if pos_z_diff_focus_tx( pos_z_diff_rx_min_index ) <= 0

    %
    image_compound(index_pos_z_rx, :) = images{ pos_z_diff_rx_min_index }( index_pos_z_rx, : );

end % for index_pos_z_rx = 1:N_pos_z_rx

end % function [ images, image_compound, image_compound_linear ] = synthesize_progressive_scan( data_RF, pos_x_focus_tx, pos_z_focus_tx, pos_z_focus_rx, positions_ctr_x, F_number_tx, F_number_rx, f_lb, f_ub, index_t0, c_ref, f_s, dyn_rx_focus )
