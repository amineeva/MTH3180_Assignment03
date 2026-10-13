%example of how to implement finite difference approximation
%for the first and second derivative of a function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%dfdx: approximation of fun'(x)
%d2fdx2: approximation of fun''(x)

function [dfdx,d2fdx2] = newton_approximate_derivative(fun,x)
    %set the step size to be tiny
    delta_x = 1e-6;
    
    %compute the function at different points near x
    f_left = fun(x-delta_x);
    f_0 = fun(x);
    f_right = fun(x+delta_x);

    %approximate the first derivative
    dfdx = (f_right-f_left)/(2*delta_x);
    
    %approximate the second derivative
    d2fdx2 = (f_right-2*f_0+f_left)/(delta_x.^2);
end