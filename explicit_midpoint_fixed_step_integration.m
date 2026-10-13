%Runs numerical integration using explicit midpoint approximation
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
function [t_list,X_list,h_avg, num_evals] = explicit_midpoint_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
t_start=tspan(1);
t_end=tspan(2);

N=ceil((t_end-t_start)/h_ref);
t_list=linspace(t_start,t_end,N);
X_list=zeros(N+1, length(X0));

t_list(1) = t_start;
X_list(1, :) = X0(:)';

num_evals=0;
h_avg = (t_end - t_start)/N;

for i = 1:N
    t_i=t_list(i);
    x_i=X_list(i,:)'; %Remember, this is a column vector
[x_new, step_evals]=explicit_midpoint_step(rate_func_in,t_i,x_i,h_avg);

num_evals=num_evals+step_evals; %Calculating total number of computes

t_list(i+1)=t_i+h_avg;
X_list(i+1,:)=x_new(:)';
end

