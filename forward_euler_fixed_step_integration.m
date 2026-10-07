%Runs numerical integration using forward Euler approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    % get the start and end times
    t_start = tspan(1);
    t_end = tspan(2);

    % find N
    N = ceil((t_end - t_start)/h_ref);

    % set vectors
    t_vector = linspace(t_start, t_end, N+1)'; % needs to be a column vector
    X_vector = zeros(N+1, length(X0));

    % initial conditions
    t_vector(1) = t_start;
    X_vector(1, :) = X0(:)'; % rows are time steps, columns are components of X

    % num evals
    num_evals = 0;
    h_avg = (t_end - t_start)/N; % same as h_ref?

    % loop for N values
    for n = 1:N
        t_n = t_vector(n);
        X_n = X_vector(n,:)'; % X_n is a column vector
        % run one forward euler step
        [XB,step_evals] = forward_euler_step(rate_func_in,t_n,X_n,h_avg);

        % update XA and time
        t_vector(n+1) = t_n + h_avg;
        X_vector(n+1,:) = XB(:)'; % XB is a column vector 

        % update num_evals
        num_evals = num_evals + step_evals;
    end

    t_list = t_vector;
    X_list = X_vector;

end
