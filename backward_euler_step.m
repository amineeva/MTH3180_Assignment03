function [XB,num_evals] = backward_euler_step(rate_func_in,t,XA,h)

    x_guess=1;
    dxtol=10^-6;
    ftol=10^-14;
    max_iter=200;
    dxmax=10^6;
    function G = backward_root_func(XA, x_guess,t, h, rate_func)
        G=XA+h*rate_func(t+h,x_guess)-x_guess;
    end
    
    backward_wrapper= @(x_guess) backward_root_func(XA, x_guess, t, h, rate_func_in);



    [XB, exit_flag, x_list] = integration_newton_solver(backward_wrapper,x_guess,dxtol,ftol,max_iter,dxmax);
    
    num_evals=2*length(x_list);
end
