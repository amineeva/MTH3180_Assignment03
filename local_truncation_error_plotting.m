%%% Local truncation error plotting %%%
% LTE -> error for a single time step
% p -> how LTE scales with the size of the time step

% clear plots
clear;
clc;


%%% Local Truncation Error: Compare Four Methods %%%

clear;
clc;
close all;

% Setup
t_ref = 0.492;
h_refs = logspace(-5, 1, 100);

% Test function 1: X(t) = cos(t)
X_ref = solution02(t_ref);

% Preallocate error arrays
lte_forward = zeros(size(h_refs));
lte_backward = zeros(size(h_refs));
lte_explicit_midpoint = zeros(size(h_refs));
lte_implicit_midpoint = zeros(size(h_refs));

% Calculate local truncation error for each step size

%%% solution 01
for i = 1:length(h_refs)

    % h = h_refs(i);
    % 
    % % Exact solution after one step
    % X_next = solution02(t_ref + h);
    % 
    % % Forward Euler
    % G_forward = X_ref + h * rate_func02(t_ref, X_ref);
    % lte_forward(i) = abs(G_forward - X_next);
    % 
    % % Explicit midpoint
    % k1 = rate_func02(t_ref, X_ref);
    % X_half = X_ref + (h/2) * k1;
    % k2 = rate_func02(t_ref + h/2, X_half);
    % G_explicit_midpoint = X_ref + h * k2;
    % lte_explicit_midpoint(i) = abs(G_explicit_midpoint - X_next);
    % 
    % % Backward Euler
    % % Solve G = X_ref + h*f(t_ref+h, G)
    % % For rate_func02(t,X) = -5*X + 5*cos(t) - sin(t),
    % % this equation can be solved directly.
    % G_backward = (X_ref + h *(5*cos(t_ref+h) - sin(t_ref+h))) / (1 + 5*h);
    % lte_backward(i) = abs(G_backward - X_next);
    % 
    % % Implicit midpoint
    % % Solve G = X_ref + h*f(t_ref+h/2, (X_ref+G)/2)
    % % For this linear rate function, solve directly.
    % t_mid = t_ref + h/2;
    % G_implicit_midpoint = (X_ref + h * (5*cos(t_mid) - sin(t_mid) - 2.5*X_ref))/ (1 + 2.5*h);
    % lte_implicit_midpoint(i) = abs(G_implicit_midpoint - X_next);

%%% solution 02
for i = 1:length(h_refs)

    h = h_refs(i);

    % Exact solution after one step
    X_next = solution02(t_ref + h);

    % System matrix for rate_func02
    A = [0, -1; 1, 0];

    % Forward Euler
    G_forward = X_ref + h * rate_func02(t_ref, X_ref);
    lte_forward(i) = norm(G_forward - X_next);

    % Explicit midpoint
    k1 = rate_func02(t_ref, X_ref);
    X_half = X_ref + (h/2) * k1;
    k2 = rate_func02(t_ref + h/2, X_half);
    G_explicit_midpoint = X_ref + h * k2;
    lte_explicit_midpoint(i) = norm(G_explicit_midpoint - X_next);

    % Backward Euler: G = X_ref + h*A*G
    G_backward = (eye(2) - h*A) \ X_ref;
    lte_backward(i) = norm(G_backward - X_next);

    % Implicit midpoint: G = X_ref + h*A*(X_ref+G)/2
    G_implicit_midpoint = ...
        (eye(2) - (h/2)*A) \ ((eye(2) + (h/2)*A)*X_ref);
    lte_implicit_midpoint(i) = norm(G_implicit_midpoint - X_next);

end


end

%%% Plot all four methods on one log-log plot %%%

figure;
hold on;

loglog(h_refs, lte_forward, 'o-','MarkerSize', 4, 'LineWidth', 1.2,'DisplayName', 'Forward Euler');
loglog(h_refs, lte_backward, 's-', 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Backward Euler');
loglog(h_refs, lte_explicit_midpoint, '^-',  'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Explicit midpoint');
loglog(h_refs, lte_implicit_midpoint, 'd-', 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Implicit midpoint');

xlabel('Step size ($h$)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Local Truncation Error ($\epsilon$)', 'Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
title('Local Truncation Error for Four Integration Methods', 'Interpreter', 'latex', 'FontSize', 18);
grid on;
box on;

%%% finding p-vals



%%% Log-log regression to estimate p %%%

% Step sizes required by the assignment
% (change the h_refs definition near the top too)
% h_refs = logspace(-5, 1, 100);

% Exclude errors too close to floating-point precision
error_floor = 1e-14;

valid_forward = lte_forward > error_floor;
valid_backward = lte_backward > error_floor;
valid_explicit_midpoint = lte_explicit_midpoint > error_floor;
valid_implicit_midpoint = lte_implicit_midpoint > error_floor;

% Fit log10(error) = p*log10(h) + intercept
fit_forward = polyfit(log10(h_refs(valid_forward)), log10(lte_forward(valid_forward)), 1);

fit_backward = polyfit(log10(h_refs(valid_backward)), log10(lte_backward(valid_backward)), 1);

fit_explicit_midpoint = polyfit(log10(h_refs(valid_explicit_midpoint)),log10(lte_explicit_midpoint(valid_explicit_midpoint)), 1);

fit_implicit_midpoint = polyfit(log10(h_refs(valid_implicit_midpoint)),log10(lte_implicit_midpoint(valid_implicit_midpoint)), 1);

% p is the slope (first coefficient)
p_forward = fit_forward(1);
p_backward = fit_backward(1);
p_explicit_midpoint = fit_explicit_midpoint(1);
p_implicit_midpoint = fit_implicit_midpoint(1);

Method = ["Forward Euler"; "Backward Euler"; ...
          "Explicit Midpoint"; "Implicit Midpoint"];

p_values = [p_forward; p_backward; ...
            p_explicit_midpoint; p_implicit_midpoint];

Results = table(Method, p_values, ...
    'VariableNames', {'Method', 'Estimated_p'});

disp(Results);



