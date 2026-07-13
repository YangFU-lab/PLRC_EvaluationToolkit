%% Procedure 4: Non-linear Fitting of Effective Solar Reflectance (ESR)
% -------------------------------------------------------------------------
% Description:
%   This script performs a non-linear fitting (ESR = a*T^4 + b*T + c) using 
%   the temperatures and known reflectances of 6 reference samples. 
%   Note: Temperature is converted to Kelvin (K) for fitting to align with 
%   the Stefan-Boltzmann law for thermal radiation (T^4).
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1. Load Input Data and Convert Units
run Procedure4_Input.m

% ▲ CRITICAL: Convert Celsius to Kelvin for physical radiation fitting (T^4)
temperature_K = temperature + 273.15; 
Time = (1:1:size(temperature_K, 1))';
n = length(Time);

% Extract reference temperatures and reflectances
T_ref = temperature_K(:, 1:6);  
R_ref = R_sol_ref;

% Extract PL-RC temperatures (Rare-earth, Quantum dots, Organic dyes)
T_PLRC_RE  = temperature_K(:, 7:9);   % [Red, Green, Blue]
T_PLRC_QD  = temperature_K(:, 10:12); % [Red, Green, Blue]
T_PLRC_dye = temperature_K(:, 13:15); % [Red, Green, Blue]

%% 2. Non-linear ESR Fitting: ESR = a*T^4 + b*T + c
disp('Performing non-linear fitting (ESR = a*T^4 + b*T + c)...');

% Reshape reference data into 1D vectors for curve fitting
R_vector = reshape(repmat(R_ref, n, 1), [], 1); 
T_vector = reshape(T_ref, [], 1); 

% Define the non-linear model
model = @(P, T) P(1) * T.^4 + P(2) * T + P(3); 
P0 = [1e-10, 0.01, 0]; % Initial guess for [a, b, c]

% Perform least squares curve fitting (Requires Optimization Toolbox)
options = optimset('Display', 'off'); % Suppress iteration output
P_fit = lsqcurvefit(model, P0, T_vector, R_vector, [], [], options);
a_nl = P_fit(1);
b_nl = P_fit(2);
c_nl = P_fit(3);

% Display fitted parameters
fprintf('--- Fitted Parameters ---\n');
fprintf('a = %e\n', a_nl);
fprintf('b = %f\n', b_nl);
fprintf('c = %f\n', c_nl);
fprintf('-------------------------\n');

%% 3. Visualize Fitted Curve vs Reference Data
figure('Name', 'Non-linear Fitting Model', 'Position', [100, 100, 600, 500]);
hold on; grid on;

% Generate fitted curve data (convert back to Celsius for plotting)
Temp_K = (20:0.01:40)' + 273.15; 
ESR_curve = model(P_fit, Temp_K);

plot(Temp_K - 273.15, ESR_curve, 'k--', 'LineWidth', 2);
plot(T_vector - 273.15, R_vector, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);

xlabel('Temperature (°C)'); ylabel('ESR');
title('Non-linear Fitted Curve and Reference Data');
legend('Fitted Curve', 'Reference Data', 'Location', 'northwest');
hold off;

%% 4. Calculate Fitted ESR for PL-RC Materials
% Calculate mean temperature for each PL-RC sample (in Kelvin)
T_ave_RE  = mean(T_PLRC_RE, 1);
T_ave_QD  = mean(T_PLRC_QD, 1);
T_ave_dye = mean(T_PLRC_dye, 1);

% Calculate corresponding ESR using the fitted model
ESR_PLRC_RE  = model(P_fit, T_ave_RE);
ESR_PLRC_QD  = model(P_fit, T_ave_QD);
ESR_PLRC_dye = model(P_fit, T_ave_dye);

%% 5. Visualization for PL-RC Materials
% Helper function for plotting to avoid redundant code
plot_PLRC_data(Temp_K, ESR_curve, T_PLRC_RE, ESR_PLRC_RE, n, 'Rare Earth-doped Phosphors');
plot_PLRC_data(Temp_K, ESR_curve, T_PLRC_QD, ESR_PLRC_QD, n, 'Quantum Dots');
plot_PLRC_data(Temp_K, ESR_curve, T_PLRC_dye, ESR_PLRC_dye, n, 'Organic Dyes');

%% --- Local Helper Function for Plotting ---
function plot_PLRC_data(Temp_K, ESR_curve, T_PLRC, ESR_PLRC, n, material_name)
    figure('Name', ['PL-RC: ', material_name], 'Position', [300, 200, 600, 500]);
    hold on; grid on;
    
    % Plot fitted curve
    plot(Temp_K - 273.15, ESR_curve, 'k--', 'LineWidth', 2);
    
    % Prepare 1D vectors for Red, Green, Blue samples
    T_R = reshape(T_PLRC(:,1), [], 1); ESR_R = reshape(repmat(ESR_PLRC(1), n, 1), [], 1);
    T_G = reshape(T_PLRC(:,2), [], 1); ESR_G = reshape(repmat(ESR_PLRC(2), n, 1), [], 1);
    T_B = reshape(T_PLRC(:,3), [], 1); ESR_B = reshape(repmat(ESR_PLRC(3), n, 1), [], 1);
    
    % Plot data points (Convert T back to Celsius)
    plot(T_R - 273.15, ESR_R, 'ro', 'MarkerSize', 4);
    plot(T_G - 273.15, ESR_G, 'go', 'MarkerSize', 4);
    plot(T_B - 273.15, ESR_B, 'bo', 'MarkerSize', 4);
    
    xlabel('Temperature (°C)'); ylabel('ESR');
    title(['Fitted Curve vs ', material_name, ' Data']);
    legend('Fitted Curve', 'Red-emissive', 'Green-emissive', 'Blue-emissive', 'Location', 'northwest');
    hold off;
end