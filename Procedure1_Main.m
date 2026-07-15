%% Procedure 1: Spectral Reconstruction of PL-RC Materials
% -------------------------------------------------------------------------
% Description:
%   This main script calculates the effective spectral reflectance & ESR
%   value of PL-RC materials by reconstructing the reflectance spectrum using 
%   filter-modified measurements and fluorescence properties.
%
% Dependencies:
%   - Procedure1_input_data.m (Must be in the same directory)
%   - Curve Fitting Toolbox (for 'fit' function)
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1. Load Input Data
run Procedure1_Input.m 

%% 2. Define Wavelength Ranges
% Define excitation and emission boundaries (Unit: nm)
lambda_ex_min = 250;  % Min wavelength of excitation
lambda_ex_max = 460;  % Max wavelength of excitation (determined by PL spectra)
lambda_em_min = 461;  % Min wavelength of emission
lambda_em_max = 800;  % Max wavelength of emission

% Generate wavelength arrays
lambda  = (250:1:2500)';                       % Full solar wavelength range
lambda1 = (lambda_ex_min:1:lambda_ex_max)';    % Excitation wavelength range
lambda2 = (lambda_em_min:1:lambda_em_max)';    % Emission wavelength range
lambda3 = (lambda_em_max + 1:1:2500)';         % Remaining wavelengths

%% 3. Data Interpolation & Visualization
disp('Interpolating spectral data...');

% 3.1 Filter-modified spectral reflectance
R_F_PL = interp1(r_filter_PL(:,1), r_filter_PL(:,2), lambda, 'linear', 'extrap');
figure('Name', 'Filter-Modified Spectral Reflectance');
plot(lambda, R_F_PL, 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('Reflectance'); title('Filter-Modified Reflectance');

% 3.2 Spectral absolute absorption @ excitation range
A_PL = interp1(Abs(:,1), Abs(:,2), lambda1, 'linear', 'extrap');
figure('Name', 'Absolute Absorption');
plot(lambda1, A_PL, 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('Absorption'); title('Absolute Absorption');

% 3.3 Spectral PLQY
PLQY = interp1(plqy(:,1), plqy(:,2), lambda1, 'linear', 'extrap');
figure('Name', 'Spectral PLQY');
plot(lambda1, PLQY, 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('PLQY'); title('Spectral PLQY');

% 3.4 Normalized emission spectrum
I_em = interp1(i_emission(:,1), i_emission(:,2), lambda2, 'linear', 'extrap');
I_em_norm = I_em / trapz(lambda2, I_em); % Area normalization

% 3.5 Spectral intensity of light source (AM 1.5)
I_sc    = interp1(i_source(:,1), i_source(:,2), lambda, 'linear', 'extrap');
I_sc_ex = interp1(i_source(:,1), i_source(:,2), lambda1, 'linear', 'extrap');
I_sc_em = interp1(i_source(:,1), i_source(:,2), lambda2, 'linear', 'extrap');

%% 4. Reconstruction of Spectral Reflectance
disp('Reconstructing effective spectral reflectance...');

% Linear fit for the filter-modified reflectance
fitr_filter_PL = fit(r_filter_PL(:,1), r_filter_PL(:,2), 'linear');

% 4.1 Corrected spectral reflectance @ excitation range
R_rec1 = fitr_filter_PL(lambda1); 

% 4.2 Corrected spectral reflectance @ emission range
% Calculate total converted photon number based on energy conservation
ConvertedPhotonNumber = trapz(lambda1, I_sc_ex .* A_PL .* PLQY .* lambda1);

R_rec2 = zeros(length(lambda2), 1);
for i = 1:length(lambda2)
    % Distribute converted photons into emission spectrum
    Temp = (ConvertedPhotonNumber * I_em_norm(i)) / (I_sc_em(i) * lambda2(i));
    R_rec2(i) = Temp + fitr_filter_PL(lambda2(i)); 
end

% 4.3 Corrected spectral reflectance @ remaining wavelengths
R_rec3 = fitr_filter_PL(lambda3); 

%% 5. Output Effective Spectral Reflectance
% Combine all ranges
R_rec = [R_rec1; R_rec2; R_rec3];

% Interpolate to standard 10 nm intervals for output
lam = (250:10:2500)';
rec = interp1(lambda, R_rec, lam, 'linear', 'extrap');

% Format output (Descending/Ascending order as needed)
OUTPUT = flipud([lam, rec]);

figure('Name', 'Effective Spectral Reflectance');
plot(lam, rec, 'k', 'LineWidth', 2);
xlabel('Wavelength (nm)'); ylabel('Effective Reflectance'); 
title('Reconstructed Effective Spectral Reflectance');

%% 6. Effective Solar Reflectance (ESR) Calculation
% Calculate the overall solar reflectance weighted by AM 1.5 spectrum
P = R_rec .* I_sc;
ESR = trapz(lambda, P) / trapz(lambda, I_sc);

fprintf('========================================\n');
fprintf('Calculated Effective Solar Reflectance (ESR): %.4f\n', ESR);
fprintf('========================================\n');
