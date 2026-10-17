
%%% Global truncation error: Four required plots %%%
clear;
clc;
close all;

%% Test Function 2
X_exact = @solution02;
X0 = solution02(0);

A = [0, -1; 1, 0]; % matrix
I = eye(2);

% %% Test Function 1
% X_exact = @solution01;
% X0 = solution01(0);
% 
% A = [-5];
% I = eye(1);

%% setup
% time interval
t0 = 0;
tf = 1;

% step sizes 
% varying N instead of h_ref case bc # of steps impacts step sizes
N_refs = unique(round(logspace(1, 4, 100))); % number of time steps
h_refs = (tf - t0) ./ N_refs; % step sizes 

num_h = length(h_refs);

%% setup arrays to store vals
gte_forward = zeros(size(h_refs));
gte_backward = zeros(size(h_refs));
gte_explicit_midpoint = zeros(size(h_refs));
gte_implicit_midpoint = zeros(size(h_refs));

calls_forward = zeros(size(h_refs));
calls_backward = zeros(size(h_refs));
calls_explicit_midpoint = zeros(size(h_refs));
calls_implicit_midpoint = zeros(size(h_refs));

%% Calculate global truncation error

for i = 1:num_h

    h = h_refs(i);
    N = N_refs(i);

    Xf = X0;
    Xb = X0;
    Xm = X0;
    Xim = X0;

    for n = 1:N

        t = t0 + (n-1)*h;

        % Forward Euler
        Xf = Xf + h*A*Xf;

        % Backward Euler
        Xb = (I - h*A) \ Xb;

        % Explicit midpoint
        k1 = A*Xm;
        Xhalf = Xm + (h/2)*k1;
        k2 = A*Xhalf;
        Xm = Xm + h*k2;

        % Implicit midpoint
        Xim = (I - (h/2)*A)\((I + (h/2)*A)*Xim);

    end

    % Global error at the final time
    Xf_exact = X_exact(tf);

    gte_forward(i) = norm(Xf - Xf_exact);
    gte_backward(i) = norm(Xb - Xf_exact);
    gte_explicit_midpoint(i) = norm(Xm - Xf_exact);
    gte_implicit_midpoint(i) = norm(Xim - Xf_exact);

    % num steps counting
    % Forward Euler: 1 evaluation per step
    % Backward Euler: 1 evaluation per step assumed
    % Explicit midpoint: 2 evaluations per step
    % Implicit midpoint: 1 evaluation per step assumed
    %
    % These are nominal counts for comparison, not literal
    % MATLAB calls in the direct linear-system implementation.
    calls_forward(i) = N;
    calls_backward(i) = N;
    calls_explicit_midpoint(i) = 2*N;
    calls_implicit_midpoint(i) = N;

end

%% fit lines (explicit)

error_floor = 1e-14;
valid_forward = gte_forward > error_floor;
valid_explicit_midpoint = gte_explicit_midpoint > error_floor;

% Global error versus step size
fit_forward_h = polyfit(log10(h_refs(valid_forward)), log10(gte_forward(valid_forward)), 1);
fit_midpoint_h = polyfit(log10(h_refs(valid_explicit_midpoint)), log10(gte_explicit_midpoint(valid_explicit_midpoint)), 1);
p_forward_h = fit_forward_h(1);
p_midpoint_h = fit_midpoint_h(1);

% Global error versus rate-function calls
fit_forward_calls = polyfit(log10(calls_forward(valid_forward)), log10(gte_forward(valid_forward)), 1);
fit_midpoint_calls = polyfit(log10(calls_explicit_midpoint(valid_explicit_midpoint)), log10(gte_explicit_midpoint(valid_explicit_midpoint)), 1);

% Positive order values: error scales as calls^(-p)
p_forward_calls = -fit_forward_calls(1);
p_midpoint_calls = -fit_midpoint_calls(1);

%% Plot 1: Global error vs. rate-function calls
% four methods

figure;
hold on;

loglog(calls_forward, gte_forward, 'o-', 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Forward Euler');

loglog(calls_backward, gte_backward, 's-','MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Backward Euler');

loglog(calls_explicit_midpoint, gte_explicit_midpoint, '^-', 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Explicit midpoint');

loglog(calls_implicit_midpoint, gte_implicit_midpoint, 'd-', 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Implicit midpoint');

xlabel('Number of Rate Function Calls', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Global Truncation Error ($\epsilon$)',  'Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
title('Global Truncation Error vs. Rate Function Calls', 'Interpreter', 'latex', 'FontSize', 18);
grid on;
box on;
hold off;

%% Plot 2: Global error vs. rate-function calls
% explicit + fit lines

calls_fit_forward = logspace(log10(min(calls_forward(valid_forward))),log10(max(calls_forward(valid_forward))), 200);

calls_fit_midpoint = logspace(log10(min(calls_explicit_midpoint(valid_explicit_midpoint))), log10(max(calls_explicit_midpoint(valid_explicit_midpoint))), 200);
forward_fit_calls = 10.^polyval(fit_forward_calls, log10(calls_fit_forward));
midpoint_fit_calls = 10.^polyval(fit_midpoint_calls, log10(calls_fit_midpoint));
figure;
hold on;
loglog(calls_forward, gte_forward, 'bo','MarkerSize', 4, 'DisplayName', 'Forward Euler Points');
loglog(calls_fit_forward, forward_fit_calls, 'b-','LineWidth', 1.5, 'DisplayName', sprintf('Forward Euler Fit: p = %.3f', p_forward_calls));
loglog(calls_explicit_midpoint, gte_explicit_midpoint, 'rs', 'MarkerSize', 4,'DisplayName', 'Explicit Midpoint Points');
loglog(calls_fit_midpoint, midpoint_fit_calls, 'r-', 'LineWidth', 1.5,  'DisplayName', sprintf('Explicit Midpoint Fit: p = %.3f', p_midpoint_calls));
xlabel('Number of Rate Function Calls', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Global Truncation Error Magnitude ($\epsilon$)','Interpreter', 'latex', 'FontSize', 18);
title('Global Truncation Error: Explicit Methods (Test Function 2)','Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
grid on;
box on;
hold off;



%% Plot 3: Global error vs. step size
% four methods

figure;
hold on;
loglog(h_refs, gte_forward, 'o-', 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Forward Euler');
loglog(h_refs, gte_backward, 's-', 'MarkerSize', 4, 'LineWidth', 1.2,'DisplayName', 'Backward Euler');
loglog(h_refs, gte_explicit_midpoint, '^-','MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Explicit midpoint');
loglog(h_refs, gte_implicit_midpoint, 'd-', 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Implicit midpoint');
xlabel('Step size ($h$)','Interpreter', 'latex', 'FontSize', 18);
ylabel('Global Truncation Error ($\epsilon$)', 'Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'best','Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
title('Global Truncation Error vs. Step Size', 'Interpreter', 'latex', 'FontSize', 18);
grid on;
box on;
hold off;

%% Plot 4: Global error vs. step size
% explicit + fit
h_fit_forward = logspace(log10(min(h_refs(valid_forward))),log10(max(h_refs(valid_forward))), 200);
h_fit_midpoint = logspace(log10(min(h_refs(valid_explicit_midpoint))), log10(max(h_refs(valid_explicit_midpoint))), 200);
forward_fit_h = 10.^polyval(fit_forward_h, log10(h_fit_forward));
midpoint_fit_h = 10.^polyval(fit_midpoint_h, log10(h_fit_midpoint));

figure;
hold on;
loglog(h_refs, gte_forward, 'bo','MarkerSize', 4,  'DisplayName', 'Forward Euler Points');
loglog(h_fit_forward, forward_fit_h, 'b-', 'LineWidth', 1.5,'DisplayName', sprintf('Forward Euler Fit: p = %.3f', p_forward_h));
loglog(h_refs, gte_explicit_midpoint, 'rs',  'MarkerSize', 4,  'DisplayName', 'Explicit Midpoint Points');
loglog(h_fit_midpoint, midpoint_fit_h, 'r-',  'LineWidth', 1.5,'DisplayName', sprintf('Explicit Midpoint Fit: p = %.3f', p_midpoint_h));
xlabel('Step size ($h$)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Global Truncation Error Magnitude ($\epsilon$)','Interpreter', 'latex', 'FontSize', 18);
title('Global Truncation Error: Explicit Methods (Test Function 2)', 'Interpreter', 'latex', 'FontSize', 18);

legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14);

set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
grid on;
box on;
hold off;


%% p values :(

error_floor = 1e-14;

% Store errors and function-call counts
errors = {gte_forward, gte_backward, gte_explicit_midpoint, gte_implicit_midpoint};

calls = {calls_forward, calls_backward, calls_explicit_midpoint, calls_implicit_midpoint};

Method = ["Forward Euler"; "Backward Euler"; "Explicit Midpoint"; "Implicit Midpoint"];

p_from_h = zeros(4,1);
p_from_calls = zeros(4,1);

for j = 1:4

    % Exclude errors near floating-point precision
    valid = errors{j} > error_floor;

    % Fit global error versus step size
    fit_h = polyfit(log10(h_refs(valid)), ...
                    log10(errors{j}(valid)), 1);

    p_from_h(j) = fit_h(1);

    % Fit global error versus function calls
    fit_calls = polyfit(log10(calls{j}(valid)), ...
                        log10(errors{j}(valid)), 1);

    % Negate slope because error decreases as calls increase
    p_from_calls(j) = -fit_calls(1);

end

Results = table(Method, p_from_h, p_from_calls, ...
    'VariableNames', {'Method', 'p_from_h', 'p_from_calls'});

disp(Results);
