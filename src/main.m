% Standard Map Analysis - Complete Executable Script
% This program analyzes the Chirikov Standard Map with functional programming
% paradigm, examining phase portraits and Lyapunov exponents

%% Define core functions
% Modulo 2π function
mod2pi = @(x) mod(x, 2*pi);

% Standard map step function
standard_map_step = @(state, K) [mod2pi(state(1) + K * sin(state(2))); mod2pi(state(2) + mod2pi(state(1) + K * sin(state(2))))];

% Function to generate a trajectory
generate_trajectory = @(initial_state, K, num_steps) trajectory_generator(initial_state, K, num_steps, standard_map_step);

% Function to generate phase portrait
generate_phase_portrait = @(K, num_trajectories, trajectory_length) phase_portrait_generator(K, num_trajectories, trajectory_length, generate_trajectory);

% Plot phase portrait
plot_phase_portrait = @(phase_portrait, K) plot_portrait(phase_portrait, K);

% Jacobian matrix (system matrix)
jacobian = @(I, theta, K) [1, K*cos(theta); 1, 1 + K*cos(theta)];

% Function to compute Lyapunov exponents
compute_lyapunov = @(K, initial_state, num_iterations, discard_iterations) lyapunov_calculator(K, initial_state, num_iterations, discard_iterations, standard_map_step, jacobian);

% Function to compute Lyapunov exponents for a range of K values
compute_lyapunov_spectrum = @(k_range, initial_state, num_iterations, discard_iterations) spectrum_calculator(k_range, initial_state, num_iterations, discard_iterations, compute_lyapunov);

%% Part (a): Phase Portraits for Different K Values
% Parameters
num_trajectories = 50;
trajectory_length = 1000;

% Random K values in specified ranges
K1 = rand * 0.6;  % K ∈ (0, 0.6]
K2 = 0.9 + rand * 0.2;  % K ∈ [0.9, 1.1]
K3 = 1.4 + rand * 0.6;  % K ∈ [1.4, 2.0]

% Generate phase portraits
phase_portrait1 = generate_phase_portrait(K1, num_trajectories, trajectory_length);
phase_portrait2 = generate_phase_portrait(K2, num_trajectories, trajectory_length);
phase_portrait3 = generate_phase_portrait(K3, num_trajectories, trajectory_length);

% Plot phase portraits
plot_phase_portrait(phase_portrait1, K1);
plot_phase_portrait(phase_portrait2, K2);
plot_phase_portrait(phase_portrait3, K3);

%% Part (b): Lyapunov Exponents
% Parameters for Lyapunov exponent calculation
num_k_values = 100;
k_range = linspace(0, 4, num_k_values);
initial_state = pi * ones(2, 1);  % Start from the middle of the phase space
num_iterations = 1000;
discard_iterations = 200;  % Discard initial transient iterations

% Compute Lyapunov spectrum
[spectrum_k, spectrum_lambda1, spectrum_lambda2] = compute_lyapunov_spectrum(k_range, initial_state, num_iterations, discard_iterations);

% Plot Lyapunov exponents
figure;
plot(spectrum_k, spectrum_lambda1, 'b-', 'LineWidth', 1.5);
hold on;
plot(spectrum_k, spectrum_lambda2, 'r-', 'LineWidth', 1.5);
plot(spectrum_k, zeros(size(spectrum_k)), 'k--');
plot(spectrum_k, spectrum_lambda1 + spectrum_lambda2, 'g-', 'LineWidth', 1);
xlabel('K');
ylabel('\lambda');
title('Lyapunov Exponents vs. K');
legend('\lambda_1', '\lambda_2', 'Zero line', '\lambda_1 + \lambda_2');
grid on;
hold off;

%% Function Definitions

% Function to generate a single trajectory
function trajectory = trajectory_generator(initial_state, K, num_steps, step_function)
trajectory = zeros(2, num_steps+1);
trajectory(:, 1) = initial_state;

for i = 1:num_steps
    trajectory(:, i+1) = step_function(trajectory(:, i), K);
end
end

% Function to generate multiple trajectories for a phase portrait
function phase_portrait = phase_portrait_generator(K, num_trajectories, trajectory_length, trajectory_generator)
phase_portrait = cell(num_trajectories, 1);

for i = 1:num_trajectories
    % Random initial conditions within [0, 2π] × [0, 2π]
    initial_state = 2 * pi * rand(2, 1);
    phase_portrait{i} = trajectory_generator(initial_state, K, trajectory_length);
end
end

% Function to plot a phase portrait
function plot_portrait(phase_portrait, K)
figure;
hold on;

for i = 1:length(phase_portrait)
    trajectory = phase_portrait{i};
    plot(trajectory(2, :), trajectory(1, :), '.', 'MarkerSize', 1);
end

xlabel('\theta');
ylabel('I');
title(['Phase Portrait for K = ', num2str(K)]);
xlim([0 2*pi]);
ylim([0 2*pi]);
grid on;
box on;
set(gca, 'XTick', 0:pi/2:2*pi);
set(gca, 'XTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
set(gca, 'YTick', 0:pi/2:2*pi);
set(gca, 'YTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
hold off;
end

% Function to calculate Lyapunov exponents for a single K value
function [lambda1, lambda2] = lyapunov_calculator(K, initial_state, num_iterations, discard_iterations, step_function, jacobian_function)
% Initialize state and Q matrix
state = initial_state;
Q = eye(2);

% Arrays to store the logarithms of the diagonal elements
log_diag1 = zeros(num_iterations, 1);
log_diag2 = zeros(num_iterations, 1);

% Discard initial transient iterations
for i = 1:discard_iterations
    state = step_function(state, K);
    Q = jacobian_function(state(1), state(2), K) * Q;
    [Q, ~] = qr(Q);
end

% Main loop for Lyapunov exponent calculation
for i = 1:num_iterations
    % Update state
    state = step_function(state, K);

    % Calculate A_n+1 = DF * Q_n
    A = jacobian_function(state(1), state(2), K) * Q;

    % QR decomposition
    [Q, R] = qr(A);

    % Store logarithms of diagonal elements
    log_diag1(i) = log(abs(R(1, 1)));
    log_diag2(i) = log(abs(R(2, 2)));
end

% Calculate Lyapunov exponents
lambda1 = sum(log_diag1) / num_iterations;
lambda2 = sum(log_diag2) / num_iterations;
end

% Function to calculate Lyapunov spectrum across a range of K values
function [spectrum_k, spectrum_lambda1, spectrum_lambda2] = spectrum_calculator(k_range, initial_state, num_iterations, discard_iterations, lyap_calculator)
spectrum_k = k_range;
num_k = length(k_range);
spectrum_lambda1 = zeros(num_k, 1);
spectrum_lambda2 = zeros(num_k, 1);

for i = 1:num_k
    K = k_range(i);
    [lambda1, lambda2] = lyap_calculator(K, initial_state, num_iterations, discard_iterations);
    spectrum_lambda1(i) = lambda1;
    spectrum_lambda2(i) = lambda2;

    % Display progress
    if mod(i, max(1, floor(num_k/10))) == 0
        fprintf('Computing Lyapunov exponents: %d%% complete\n', round(100*i/num_k));
    end
end
end
