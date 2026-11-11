%% Procedure 3: Linear Fitting of Effective Solar Reflectance (ESR)
% -------------------------------------------------------------------------
% Description:
%   This script fits a linear relationship (ESR = a*T + b) using the
%   temperatures and known solar reflectances of six reference samples.
%   It then estimates ESR for three PL-RC materials.
%
%   Goodness-of-fit metrics are calculated ONLY for the selected
%   steady-state period:
%       1. R-squared
%       2. RMSE = sqrt(SSE / n_ref)
%       3. Maximum absolute residual
%
% Dependencies:
%   - Procedure3_Input.m (must be in the same directory)
%
% Notes:
%   - Temperatures are expressed in degrees Celsius.
%   - Reflectances and ESR values use a fractional scale, not percentages.
%   - The steady-state interval must be selected from the measured data.
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1. Load Input Data
run Procedure3_Input.m

% Validate input dimensions
assert(isnumeric(temperature) && size(temperature, 2) >= 9, ...
    'temperature must contain at least nine columns.');

assert(all(all(isfinite(temperature(:, 1:9)))), ...
    'The first nine temperature columns must contain only finite values.');

% Measurement index; this is not elapsed time unless a sampling interval
% is explicitly provided.
Time = (1:size(temperature, 1))';

% Reference samples (columns 1 to 6)
T_ref = temperature(:, 1:6);
R_ref = reshape(R_sol_ref, 1, []);

assert(numel(R_ref) == 6, ...
    'R_sol_ref must contain six reference reflectances.');

assert(all(isfinite(R_ref)), ...
    'Reference reflectances must contain only finite values.');

assert(sum((R_ref - mean(R_ref)).^2) > 0, ...
    'Reference reflectances must not all be identical.');

% PL-RC samples (columns 7 to 9)
T_PLRC_B_RE = temperature(:, 7);
T_PLRC_G_RE = temperature(:, 8);
T_PLRC_R_RE = temperature(:, 9);

%% 2. Continuous Linear ESR Fitting: ESR = a*T + b
disp('Performing continuous linear fitting over time...');

% Define temperature range for plotting the fitted lines
Temp = (10:0.1:22)';

% Define the analysis range (default: all measurements)
n1 = 1;
n2 = length(Time);
num_points = n2 - n1 + 1;

% Preallocate ESR arrays
ESR_b = zeros(num_points, 1);
ESR_g = zeros(num_points, 1);
ESR_r = zeros(num_points, 1);

figure('Name', 'Real-time Linear Fitting', ...
    'Position', [100, 100, 600, 500]);
hold on;
grid on;

for i = n1:n2
    idx = i - n1 + 1;

    assert(max(T_ref(i, :)) > min(T_ref(i, :)), ...
        'Reference temperatures are identical at measurement %d.', i);

    % Fit the six reference samples at the current time point
    p = polyfit(T_ref(i, :), R_ref, 1);

    % Estimate ESR for the three PL-RC materials
    ESR_b(idx) = polyval(p, T_PLRC_B_RE(i));
    ESR_g(idx) = polyval(p, T_PLRC_G_RE(i));
    ESR_r(idx) = polyval(p, T_PLRC_R_RE(i));

    % Plot the fitted line and reference data
    plot(Temp, polyval(p, Temp), ...
        'Color', [0.7 0.7 0.7]);

    plot(T_ref(i, :), R_ref, ...
        'ko', 'MarkerSize', 4);
end

xlabel('Temperature (°C)');
ylabel('ESR');
title('Real-time ESR vs Temperature Fitting');
hold off;

%% 3. Visualization of Fitted ESR over Time
figure('Name', 'Fitted ESR of PL-RC Materials', ...
    'Position', [750, 100, 600, 500]);

plot(Time(n1:n2), ESR_b, 'b', 'LineWidth', 1.5);
hold on;
plot(Time(n1:n2), ESR_g, 'g', 'LineWidth', 1.5);
plot(Time(n1:n2), ESR_r, 'r', 'LineWidth', 1.5);

xlabel('Measurement index');
ylabel('Fitted ESR');
title('Fitted ESR Dynamics of PL-RC Materials');

legend('Blue-emissive', 'Green-emissive', 'Red-emissive', ...
    'Location', 'best');

grid on;
hold off;

% Statistics over the full analysis range selected in Section 2
ESR_mean = [mean(ESR_b), mean(ESR_g), mean(ESR_r)];
ESR_min  = [min(ESR_b), min(ESR_g), min(ESR_r)];
ESR_max  = [max(ESR_b), max(ESR_g), max(ESR_r)];

fprintf('\n--- Statistical Results over the Analysis Range ---\n');

fprintf('Blue  - Mean: %.4f, Min: %.4f, Max: %.4f\n', ...
    ESR_mean(1), ESR_min(1), ESR_max(1));

fprintf('Green - Mean: %.4f, Min: %.4f, Max: %.4f\n', ...
    ESR_mean(2), ESR_min(2), ESR_max(2));

fprintf('Red   - Mean: %.4f, Min: %.4f, Max: %.4f\n', ...
    ESR_mean(3), ESR_min(3), ESR_max(3));

fprintf('---------------------------------------------------\n');

%% 4. Linear Fitting Within a Selected Steady-State Period
disp('Performing linear fitting for the selected steady-state period...');

% IMPORTANT:
% Adjust m1 and m2 according to the actual steady-state interval.
% These values are row indices, not elapsed times.
m1 = 586;
m2 = 595;

assert(isscalar(m1) && isscalar(m2) && ...
    isfinite(m1) && isfinite(m2) && ...
    m1 == fix(m1) && m2 == fix(m2) && ...
    m1 >= 1 && m2 > m1 && m2 <= size(temperature, 1), ...
    'm1 and m2 must define a valid interval containing at least two rows.');

% Mean temperature of each reference sample during the selected period
T_ave = mean(T_ref(m1:m2, :), 1);

assert(max(T_ave) > min(T_ave), ...
    'Mean reference temperatures must not all be identical.');

% Fit the six mean reference temperatures
p_ave = polyfit(T_ave, R_ref, 1);

% Display the fitted equation
equation_str = sprintf( ...
    'ESR = %.6f * T %+.6f', p_ave(1), p_ave(2));

fprintf('Fitted equation: %s\n', equation_str);

%% 4.1. Goodness-of-Fit for the Steady-State Calibration
% All metrics below are evaluated using the SIX mean reference points.

% Predicted reflectance at each mean reference temperature
R_fit_ave = polyval(p_ave, T_ave);

% Signed residual: measured reflectance minus fitted reflectance
residual_ave = R_ref - R_fit_ave;

% Sum of squared residuals and total sum of squares
n_ref = numel(R_ref);
SSE_ave = sum(residual_ave.^2);
SST_ave = sum((R_ref - mean(R_ref)).^2);

% Goodness-of-fit metrics
R2_ave = 1 - SSE_ave / SST_ave;
RMSE_ave = sqrt(SSE_ave / n_ref);
MaxResidual_ave = max(abs(residual_ave));

% Display the overall metrics
fprintf('R-squared: %.5f\n', R2_ave);
fprintf('RMSE: %.5f\n', RMSE_ave);
fprintf('Maximum absolute residual: %.5f\n', MaxResidual_ave);
fprintf('--------------------------------\n');

%% 4.2. Visualization of the Steady-State Calibration
% Individual time-point measurements within the selected interval
R_datapoint = reshape( ...
    repmat(R_ref, m2 - m1 + 1, 1), [], 1);

T_datapoint = reshape(T_ref(m1:m2, :), [], 1);

figure('Name', 'Steady-State Fitting', ...
    'Position', [400, 300, 700, 500]);
hold on;
grid on;

% Fitted calibration line
hFit = plot(Temp, polyval(p_ave, Temp), ...
    'r--', 'LineWidth', 2);

% Individual measurements within the selected interval
hRaw = plot(T_datapoint, R_datapoint, ...
    'bo', 'MarkerSize', 4);

% Six mean reference points used for fitting
hMean = plot(T_ave, R_ref, ...
    'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 7);

xlabel('Temperature (°C)');
ylabel('ESR');

title({ ...
    'Fitted Curve vs Original Data (Stable Period)', ...
    sprintf('R^2 = %.4f; RMSE = %.4f; Max. |residual| = %.4f', ...
    R2_ave, RMSE_ave, MaxResidual_ave)});

legend([hFit, hRaw, hMean], ...
    {'Fit to mean temperatures', ...
     'Individual time-point measurements', ...
     'Mean reference temperatures'}, ...
    'Location', 'best');

% Display the fitted equation in the figure
text(0.05, 0.08, equation_str, ...
    'Units', 'normalized', ...
    'FontSize', 11, ...
    'Color', 'r', ...
    'FontWeight', 'bold', ...
    'Interpreter', 'none');

hold off;