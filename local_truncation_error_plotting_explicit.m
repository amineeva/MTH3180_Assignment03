
%%% Local Truncation Error Experiment %%%

clear;
clc;

% Setup
t_ref = 0.492;
h_refs = logspace(-5, 1, 100);

% Test functions
rate_func = @rate_func01; %switch to @rate_func02 for second test
solution = @solution01; %switch to @solution02 for second test
label = 'Test Function 1';

% Columns: Forward Euler, Explicit Midpoint, Exact Solution Change
p_results = zeros(2, 3);

% Exact solution at the reference time
X_ref = solution(t_ref);

% Preallocate arrays
lte_euler = zeros(size(h_refs));
lte_midpoint = zeros(size(h_refs));
delta_X = zeros(size(h_refs));

% Calculate errors for each step size
for i = 1:length(h_refs)

    h = h_refs(i);

    % Exact solution at the end of one step
    X_next = solution(t_ref + h);

    % Forward Euler
    G_euler = X_ref + h * rate_func(t_ref, X_ref);

    % Explicit midpoint
    k1 = rate_func(t_ref, X_ref);
    X_half = X_ref + (h/2) * k1;
    k2 = rate_func(t_ref + h/2, X_half);
    G_midpoint = X_ref + h * k2;

    % Norms work for both scalar and vector solutions
    lte_euler(i) = norm(G_euler - X_next);
    lte_midpoint(i) = norm(G_midpoint - X_next);
    delta_X(i) = norm(X_next - X_ref);

end

% % Exclude errors near floating-point precision from regression.
% % All data points are still retained for plotting.
error_floor = 1e-14;

valid_euler = lte_euler > error_floor;
valid_midpoint = lte_midpoint > error_floor;
valid_delta = delta_X > error_floor;

% valid_euler = lte_euler;
% valid_midpoint = lte_midpoint;
% valid_delta = delta_X;

% Log-log regression
% log10(error) = p*log10(h) + intercept
fit_euler = polyfit(log10(h_refs(valid_euler)), log10(lte_euler(valid_euler)), 1);

fit_midpoint = polyfit(log10(h_refs(valid_midpoint)), log10(lte_midpoint(valid_midpoint)), 1);

fit_delta = polyfit(log10(h_refs(valid_delta)),log10(delta_X(valid_delta)), 1);

% Save the fitted slopes
p_results(1, :) = [fit_euler(1), fit_midpoint(1), fit_delta(1)];

% Values for smooth fit lines
h_fit = logspace(-5, 1, 200);

euler_fit = 10.^polyval(fit_euler, log10(h_fit));
midpoint_fit = 10.^polyval(fit_midpoint, log10(h_fit));
delta_fit = 10.^polyval(fit_delta, log10(h_fit));

% Plot data and regression lines
figure;
hold on;

loglog(h_refs, lte_euler, 'bo','MarkerSize', 4, 'DisplayName', 'Forward Euler Points');
loglog(h_fit, euler_fit, 'b-','LineWidth', 1.5,'DisplayName', sprintf('Forward Euler Fit: p = %.3f', fit_euler(1)));
loglog(h_refs, lte_midpoint, 'rs','MarkerSize', 4, 'DisplayName', 'Explicit Midpoint Points');
loglog(h_fit, midpoint_fit, 'r-', 'LineWidth', 1.5,'DisplayName', sprintf('Explicit Midpoint Fit: p = %.3f', fit_midpoint(1)));
loglog(h_refs, delta_X, 'k^','MarkerSize', 4, 'DisplayName', '|X(t+h) - X(t)| data');
loglog(h_fit, delta_fit, 'k--','LineWidth', 1.5, 'DisplayName', sprintf('Exact change Fit: p = %.3f', fit_delta(1)));
xlabel('Step size ($h$)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Error magnitude');
ylabel('Local Truncation Error Magnitude ($\epsilon$)', 'Interpreter', 'latex', 'FontSize', 18);
title('Local Truncation Error: Explicit Methods (Test Func 1)', 'Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');

grid on;
box on;
hold off;


% Results table
Results = table( 'labels',  p_results(:,1), p_results(:,2), p_results(:,3), 'VariableNames', {'TestFunction','ForwardEuler_p', 'ExplicitMidpoint_p', 'ExactChange_p'});

disp(Results);


