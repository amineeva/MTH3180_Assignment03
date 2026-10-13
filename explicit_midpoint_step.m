function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    %your code here
    % derivative at current point for tangent line
    mid_t=t+h/2;
    mid_x=XA+.5*h*rate_func_in(t,XA);

    XB=XA+h*rate_func_in(mid_t,mid_x);

    num_evals = 2;
    
end
