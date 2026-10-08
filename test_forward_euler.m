% test forward_euler_fixed_step_integration

% %%% Testing solution01
% tspan = [0, 10];
% h_refs = [0.55, 0.5, 0.05];
% X0 = 1;
% 
% figure;
% hold on;
% 
% % exact solution
% t_exact = linspace(tspan(1), tspan(2), 1000);
% plot(t_exact, solution01(t_exact), 'k', 'LineWidth', 2.5, 'DisplayName', 'Exact');
% 
% % using Euler solver
% for h_ref = h_refs
%     [t_list, X_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
%     plot(t_list, X_list, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', sprintf('$h_{\\mathrm{ref}} = %.2f$', h_ref));
% end
% 
% % plotting
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('$X(t)$', 'Interpreter', 'latex', 'FontSize', 18);
% legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14);
% set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
% title('Forward Euler Approximation of $X(t)=\cos(t)$', 'Interpreter', 'latex', 'FontSize', 18);
% grid on;
% box on;
% 
% hold off;


%%% Testing solution02
tspan = [0, 10];
h_refs = [0.55, 0.5, 0.05];
X0 = [1; 0];

figure;
hold on;

% exact solution
t_exact = linspace(tspan(1), tspan(2), 1000);
X_exact = solution02(t_exact);
plot(t_exact, X_exact(1,:), 'k', 'LineWidth', 2.5, 'DisplayName', 'Exact: $X_1(t)=\cos(t)$');
plot(t_exact, X_exact(2,:), 'k--', 'LineWidth', 2.5, 'DisplayName', 'Exact: $X_2(t)=\sin(t)$');

% using Euler solver
for h_ref = h_refs
    [t_list, X_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func02, tspan, X0, h_ref);
    plot(t_list, X_list(1,:), 'o-', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', sprintf('$X_1$, $h_{\\mathrm{ref}}=%.2f$', h_ref));
    plot(t_list, X_list(2,:), 's-', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', sprintf('$X_2$, $h_{\\mathrm{ref}}=%.2f$', h_ref));
end

% plotting
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$X(t)$', 'Interpreter', 'latex', 'FontSize', 18);
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
title('Forward Euler Approximation of $\mathbf{X}(t)$', 'Interpreter', 'latex', 'FontSize', 18);
grid on;
box on;

% hold off;


%%% Example 1
function dXdt = rate_func01(t, X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end


%%% Example 2
function dXdt = rate_func02(t, X)
    dXdt = [0,-1;1,0]*X;
end

function X = solution02(t)
    X = [cos(t);sin(t)];
end
