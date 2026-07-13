%% Procedure 2: Spectral Combination for Effective Spectral Reflectance
% -------------------------------------------------------------------------
% Description:
%   This main script calculates the effective spectral reflectance (ESR) 
%   by combining spectral reflectance data obtained from different light 
%   sources (DH lamp, sun simulator, and conventional spectrometer) across 
%   their respective valid wavelength ranges.
%
% Dependencies:
%   - Procedure2_input_data.m (Must be in the same directory)
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1. Load Input Data
run Procedure2_Input.m 

%% 2. Define Wavelength Ranges
% Define excitation, emission, and remaining boundaries (Unit: nm)
lambda_ex_min = 250;  % Min wavelength of excitation
lambda_ex_max = 450;  % Max wavelength of excitation
lambda_em_min = 451;  % Min wavelength of emission (determined by PL spectra)
lambda_em_max = 800;  % Max wavelength of emission

% Generate wavelength arrays for different segments
lambda1 = (lambda_ex_min:1:lambda_ex_max)';    % Excitation range
lambda2 = (lambda_em_min:1:lambda_em_max)';    % Emission range
lambda3 = (lambda_em_max + 1:1:2500)';         % Remaining wavelengths

% Full wavelength range
lambda = [lambda1; lambda2; lambda3];

%% 3. Data Interpolation
disp('Interpolating spectral data for combination...');

% Interpolate data to match the defined wavelength ranges
R1_DHlamp = interp1(r_DHlamp(:,1), r_DHlamp(:,2), lambda1, 'linear', 'extrap');
R2_sunsim = interp1(r_sunsim(:,1), r_sunsim(:,2), lambda2, 'linear', 'extrap');
R3_pseudo = interp1(r_pseudo(:,1), r_pseudo(:,2), lambda3, 'linear', 'extrap');

%% 4. Effective Spectral Reflectance Combination
% Combine the three segments into a single effective reflectance spectrum
R_comb = [R1_DHlamp; R2_sunsim; R3_pseudo];

% Visualization of the combined spectrum vs original data
figure('Name', 'Spectral Combination Overview', 'Position', [100, 100, 800, 500]);
semilogy(r_DHlamp(:,1), r_DHlamp(:,2), 'b--', 'LineWidth', 1.5); hold on;
semilogy(r_sunsim(:,1), r_sunsim(:,2), 'y:', 'LineWidth', 2);
semilogy(r_pseudo(:,1), r_pseudo(:,2), 'k-', 'LineWidth', 1);
semilogy(lambda, R_comb, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');

% Formatting the plot
xlabel('Wavelength (nm)'); 
ylabel('Spectral Reflectance (Log Scale)'); 
title('Combination of Spectral Reflectance');
legend('DH Lamp Data', 'Sun Simulator Data', 'Conventional Spectrometer Data', ...
       'Combined Effective Reflectance', 'Location', 'best');
grid on;

%% 5. Effective Solar Reflectance (ESR) Calculation
disp('Calculating Effective Solar Reflectance...');

% Interpolate the solar source spectrum to match the combined wavelength array
I_sc = interp1(i_source(:,1), i_source(:,2), lambda, 'linear', 'extrap');

% Calculate the overall solar reflectance weighted by AM 1.5 spectrum
P = R_comb .* I_sc;
ESR = trapz(lambda, P) / trapz(lambda, I_sc);

% Output the final result
fprintf('========================================\n');
fprintf('Calculated Effective Solar Reflectance (ESR): %.4f\n', ESR);
fprintf('========================================\n');