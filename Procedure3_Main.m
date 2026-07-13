%% Procedure 3: Linear Fitting of Effective Solar Reflectance (ESR)
% -------------------------------------------------------------------------
% Description:
%   This script performs a linear fitting (ESR = a*T + b) using the 
%   temperatures and known reflectances of 6 reference samples. It then 
%   estimates the ESR for three PL-RC materials (Blue, Green, Red emissive) 
%   based on their measured temperatures.
%
% Dependencies:
%   - Procedure3_input_data.m (Must be in the same directory)
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1. Load Input Data
run Procedure3_Input.m

% Extract time array and temperature columns
Time = (1:1:size(temperature, 1))';

% Reference samples (1 to 6)
T_ref = temperature(:, 1:6); 
R_ref = R_sol_ref;

% PL-RC samples (Rare-earth doped phosphors: Blue, Green, Red)
T_PLRC_B_RE = temperature(:, 7);
T_PLRC_G_RE = temperature(:, 8);
T_PLRC_R_RE = temperature(:, 9);

%% 2. Continuous Linear ESR Fitting: ESR = a*T + b
disp('Performing continuous linear fitting over time...');

% Define temperature range for plotting the fitted lines
Temp = (10:0.1:22)';

% Define the time range for analysis (default: all time points)
n1 = 1;
n2 = length(Time);
num_points = n2 - n1 + 1;

% Preallocate arrays for speed and memory efficiency
ESR_b = zeros(num_points, 1);
ESR_g = zeros(num_points, 1);
ESR_r = zeros(num_points, 1);

figure('Name', 'Real-time Linear Fitting', 'Position', [100, 100, 600, 500]);
hold on; grid on;
for i = n1:n2
    idx = i - n1 + 1;
    
    % Linear fitting for the current timepoint
    p = polyfit(T_ref(i, :), R_ref, 1); 
    
    % Calculate fitted ESR for the three PL-RC materials
    ESR_b(idx) = polyval(p, T_PLRC_B_RE(i)); 
    ESR_g(idx) = polyval(p, T_PLRC_G_RE(i)); 
    ESR_r(idx) = polyval(p, T_PLRC_R_RE(i)); 
    
    % Plot the fitting line and original reference data points
    plot(Temp, polyval(p, Temp), 'Color', [0.7 0.7 0.7 0.5]); % Gray transparent lines
    plot(T_ref(i, :), R_ref, 'ko', 'MarkerSize', 4);
end
xlabel('Temperature (°C)'); ylabel('ESR');
title('Real-time ESR vs Temperature Fitting');
hold off;

%% 3. Visualization of Fitted ESR over Time
figure('Name', 'Fitted ESR of PL-RC Materials', 'Position', [750, 100, 600, 500]);
plot(Time(n1:n2), ESR_b, 'b', 'LineWidth', 1.5); hold on;
plot(Time(n1:n2), ESR_g, 'g', 'LineWidth', 1.5);
plot(Time(n1:n2), ESR_r, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Fitted ESR');
title('Fitted ESR Dynamics of PL-RC Materials');
legend('Blue-emissive', 'Green-emissive', 'Red-emissive', 'Location', 'best');
grid on;

% Calculate statistical values
ESR_mean = [mean(ESR_b), mean(ESR_g), mean(ESR_r)];
ESR_min  = [min(ESR_b), min(ESR_g), min(ESR_r)];
ESR_max  = [max(ESR_b), max(ESR_g), max(ESR_r)];

fprintf('--- Statistical Results of Fitted ESR ---\n');
fprintf('Blue  - Mean: %.4f, Min: %.4f, Max: %.4f\n', ESR_mean(1), ESR_min(1), ESR_max(1));
fprintf('Green - Mean: %.4f, Min: %.4f, Max: %.4f\n', ESR_mean(2), ESR_min(2), ESR_max(2));
fprintf('Red   - Mean: %.4f, Min: %.4f, Max: %.4f\n', ESR_mean(3), ESR_min(3), ESR_max(3));
fprintf('-----------------------------------------\n');

%% 4. Linear Fitting Within a Selected Steady-State Period
disp('Performing linear fitting for the selected steady-state period...');

% Define the specific time window for steady-state analysis 
% ▲ CAUTION: Users must adjust m1 and m2 based on their actual steady-state data!
m1 = 586; 
m2 = 595; 

% Average temperature of reference samples during the selected period
T_ave = mean(T_ref(m1:m2, :), 1);
p_ave = polyfit(T_ave, R_ref, 1);

% Display fitted equation
equation_str = sprintf('Fitted Equation: ESR = %.4f * T + %.4f', p_ave(1), p_ave(2));
disp(equation_str);

% Prepare data points for visualization
R_datapoint = reshape(repmat(R_ref, m2 - m1 + 1, 1), [], 1); 
T_datapoint = reshape(T_ref(m1:m2, :), [], 1); 

figure('Name', 'Steady-State Fitting', 'Position', [400, 300, 600, 500]);
hold on; grid on;
% Plot fitted curve
plot(Temp, polyval(p_ave, Temp), 'r--', 'LineWidth', 2);
% Plot original data points
plot(T_datapoint, R_datapoint, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); 

xlabel('Temperature (°C)'); ylabel('ESR');
title('Fitted Curve vs Original Data (Steady-State Period)');
legend('Average Fitted Curve', 'Original Data Points', 'Location', 'northwest');
% Add equation text to the plot
text(min(Temp)+1, max(R_ref)-0.05, equation_str, 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
hold off;