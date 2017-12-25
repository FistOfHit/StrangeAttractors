function [lyap_exponents] = Lyapunov_exponents(a, b, c, force1, force2, irr_freq)
% Find all lyapunov exponents of duffing oscillator system with or 10 dimensional method.
%
% Inputs:
% a = -1*Damping strength
% b = Linear stiffness
% c = non-linearity in restoring force
% force1 = amplitude of force in rational frequency driving force (2*pi Hz)
% force2 = amplitude of force in irrational frequecy driving force (2*pi*irr_freq Hz)
% irr_freq = value of irrational driving frequency
%
% Outputs:
% lyap_exponents = array of all 4 estimations of lyapunov exponents of system
%
% This obscenely long and complicated method of finding the lyapunov exponents is explained much
% better in the final report available in the repository, but as a summary: We clone the system
% into 5 systems and solve simultaneously, with a perturbation in a different dimenison in each
% system alongside a regular system. We then perform a slightly modified version of Gram-Schimdt to
% ensure the flows are seperate and hence finding each individual exponent. Its reccomended not to
% look too much at the blocks of code to avoid discomfort and trouble with eye sight. Here, the
% first skip_iter iterations are skipped to absolutelty ensure we are only working with the long
% term "stabilised" solution that isn't unwantedly affected by initial conditios (which are default
% set to 1, 0 here, as the lyapunov exponent is unaffected by them).


    % Initialising values needed
    step_size = 100;
    t_step= (2*pi)/step_size;
    tau_step = t_step*irr_freq;
    
    % Initial conditons for orbits and tangent vectors
    p = 1;
    s = 0;
    z = 0;
    z_tau = 0;
    
    % System perturbed in x
    dp1 = 10^(-2);
    ds1 = 0;
    dz1 = 0;
    dz_tau1 = 0;
    
    % System perturbed in dx/dt
    dp2 = 0;
    ds2 = 10^(-2);
    dz2 = 0;
    dz_tau2 = 0;
    
    % System perturbed in time
    dp3 = 0;
    ds3 = 0;
    dz3 = 10^(-2);
    dz_tau3 = 0;
    
    % System perturbed in tau
    dp4 = 0;
    dv4 = 0;
    dz4 = 0;
    dz_tau4 = 10^(-2);
    
    % Error size we will keep our tangent line distance at
    err_size1 = norm([dp1, ds1, dz1, dz_tau1]);
    err_size2 = norm([dp2, ds2, dz2, dz_tau2]);
    err_siez3 = norm([dp3, ds3, dz3, dz_tau3]);
    err_siez4 = norm([dp4, dv4, dz4, dz_tau4]);
    
    % Counters to track progression through driving cycle
    t_step_count    = 0;
    tau_step_count = 0;
    
    skip_iter = 10^4;
    num_iter = 10^5;
    % Iterations to skip transients
    for i = 1 : skip_iter 
        
        % Each of the 4 steps of RK4 on each of 5 systems
        %---------------------------------------------------------------
        k1 = s;
        l1 = a*s + b*p + c*p^3 + force1*cos(z) + force2*cos(z_tau);
  
        m1 = ds1;
        n1 = a*ds1 + b*dp1 + 3*c*p^2*dp1 - force1*sin(z)*dz1 - force2*sin(z_tau)*dz_tau1;
        
        f1 = ds2;
        g1 = a*ds2 + b*dp2 + 3*c*p^2*dp2 - force1*sin(z)*dz2 - force2*sin(z_tau)*dz_tau2;
        
        h1 = ds3;
        j1 = a*ds3 + b*dp3 + 3*c*p^2*dp3 - force1*sin(z)*dz3 - force2*sin(z_tau)*dz_tau3;
        
        u1 = dv4;
        v1 = a*dv4 + b*dp4 + 3*c*p^2*dp4 - force1*sin(z)*dz4 - force2*sin(z_tau)*dz_tau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        k2 = s + 0.5*t_step*l1;
        l2 = a*(s + 0.5*t_step*l1) + b*(p + 0.5*t_step*k1) + c*(p + 0.5*t_step*k1)^3 +...
             force1*cos(z + 0.5*t_step) + force2*cos(z_tau + 0.5*tau_step);
        
        m2 = ds1 + 0.5*t_step*n1;
        n2 = a*(ds1 + 0.5*t_step*n1) + b*(dp1 + 0.5*t_step*m1) +...
             3*c*(p + 0.5*t_step*k1)^2 * (dp1 + 0.5*t_step*m1) -...
             force1*sin(z + 0.5*t_step)*dz1 - force2*sin(z_tau + 0.5*tau_step)*dz_tau1;
        
        f2 = ds2 + 0.5*t_step*g1;
        g2 = a*(ds2 + 0.5*t_step*g1) + b*(dp2 + 0.5*t_step*f1) +...
             3*c*(p + (0.5*t_step*k1))^2 * (dp2 + 0.5*t_step*f1) -...
             force1*sin(z + 0.5*t_step)*dz2 - force2*sin(z_tau + 0.5*tau_step)*dz_tau2;
        
        h2 = ds3 + (0.5*t_step*j1);
        j2 = a*(ds3 + 0.5*t_step*j1) + b*(dp3 + 0.5*t_step*h1) +...
             3*c*(p + 0.5*t_step*k1)*(p + 0.5*t_step*k1)*(dp3 + 0.5*t_step*h1) -...
             force1*sin(z + 0.5*t_step)*dz3 - force2*sin(z_tau + 0.5*tau_step)*dz_tau3;
        
        u2 = dv4 + 0.5*t_step*v1;
        v2 = a*(dv4 + 0.5*t_step*v1) + b*(dp4 + 0.5*t_step*u1) +...
             3*c*(p + 0.5*t_step*k1)^2 * (dp4 + 0.5*t_step*u1) -...
             force1*sin(z + 0.5*t_step)*dz4 - force2*sin(z_tau + 0.5*tau_step)*dz_tau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        k3 = s + 0.5*t_step*l2;
        l3 = a*(s + 0.5*t_step*l2) + b*(p + 0.5*t_step*k2) + c*(p + 0.5*t_step*k2)^3 +...
             force1*cos(z + 0.5*t_step) + force2*cos(z_tau + 0.5*tau_step);
        
        m3 = ds1 + 0.5*t_step*n2;
        n3 = a*(ds1 + 0.5*t_step*n2) + b*(dp1 + 0.5*t_step*m2) +...
             3*c*(p + 0.5*t_step*k2)^2 * (dp1 + 0.5*t_step*m2) -...
             force1*sin(z + 0.5*t_step)*dz1 - force2*sin(z_tau + 0.5*tau_step)*dz_tau1;
        
        f3 = ds2 + 0.5*t_step*g2;
        g3 = a*(ds2 + 0.5*t_step*g2) + b*(dp2 + 0.5*t_step*f2) +...
             3*c*(p + 0.5*t_step*k2)*(p + 0.5*t_step*k2)*(dp2 + 0.5*t_step*f2) -...
             force1*sin(z + 0.5*t_step)*dz2 - force2*sin(z_tau + 0.5*tau_step)*dz_tau2;
        
        h3 = ds3 + 0.5*t_step*j2;
        j3 = a*(ds3 + 0.5*t_step*j2) + b*(dp3 + 0.5*t_step*h2) +...
            3*c*(p + 0.5*t_step*k2)*(p + 0.5*t_step*k2)*(dp3 + 0.5*t_step*h2) -...
            force1*sin(z + 0.5*t_step)*dz3 - force2*sin(z_tau + 0.5*tau_step)*dz_tau3;
        
        u3 = dv4 + 0.5*t_step*v2;
        v3 = a*(dv4 + 0.5*t_step*v2) + b*(dp4 + 0.5*t_step*u2) +...
             3*c*(p + 0.5*t_step*k2)*(p + 0.5*t_step*k2)*(dp4 + 0.5*t_step*u2) -...
             force1*sin(z + 0.5*t_step)*dz4 - force2*sin(z_tau + 0.5*tau_step)*dz_tau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        k4 = s + t_step*l3;       
        l4 = a*(s + t_step*l3) + b*(p + t_step*k3) + c*(p + t_step*k3)^3 +...
             force1*cos(z + t_step) + force2*cos(z_tau + tau_step);
       
        m4 = ds1 + t_step*n3;
        n4 = a*(ds1 + t_step*n3) + b*(dp1 + t_step*m3) +...
            3*c*(p + t_step*k3)^2 * (dp1 + t_step*m3) -...
            force1*sin(z + t_step)*dz1 - force2*sin(z_tau + tau_step)*dz_tau1;
        
        f4 = ds2 + t_step*g3;
        g4 = a*(ds2 + t_step*g3) + b*(dp2 + t_step*f3) +...
             3*c*(p + t_step*k3)^2 * (dp2 + (t_step*f3)) -...
             force1*sin(z + t_step)*dz2 - force2*sin(z_tau + tau_step)*dz_tau2;
        
        h4 = ds3 + t_step*j3;
        j4 = a*(ds3 + t_step*j3) + b*(dp3 + t_step*h3) +...
            3*c*(p + t_step*k3)^2 * (dp3 + t_step*h3) -...
            force1*sin(z + t_step)*dz3 - force2*sin(z_tau + tau_step)*dz_tau3;
        
        u4 = dv4 + t_step*v3;
        v4 = a*(dv4 + t_step*v3) + b*(dp4 + t_step*u3) +...
             3*c*(p + t_step*k3)*(p + t_step*k3)*(dp4 + t_step*u3) -...
             force1*sin(z + t_step)*dz4 - force2*sin(z_tau + tau_step)*dz_tau4;
        %---------------------------------------------------------------
        
        p = p + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*t_step;
        s = s + (1/6)*(l1 + 2*l2 + 2*l3 + l4)*t_step;
        
        dp1 = dp1 + (1/6)*(m1 + 2*m2 + 2*m3 + m4)*t_step;
        ds1 = ds1 + (1/6)*(n1 + 2*n2 + 2*n3 + n4)*t_step;

        dp2 = dp2 + (1/6)*(force1 + 2*f2 + 2*f3 + f4)*t_step;
        ds2 = ds2 + (1/6)*(g1 + 2*g2 + 2*g3 + g4)*t_step;

        dp3 = dp3 + (1/6)*(h1 + 2*h2 + 2*h3 + h4)*t_step;
        ds3 = ds3 + (1/6)*(J1 + 2*j2 + 2*j3 + j4)*t_step;
        
        dp4 = dp4 + (1/6)*(u1 + 2*u2 + 2*u3 + u4)*t_step;
        dv4 = dv4 + (1/6)*(v1 + 2*v2 + 2*v3 + v4)*t_step;
        
        z = z + t_step;
        t_step_count = t_step_count + 1;
        
        z_tau = z_tau + tau_step;
        tau_step_count = tau_step_count + 1;
        
        if t_step_count == step_size
            z = 0;
            t_step_count = 0;
        end
        
        if tau_step_count == step_size
            z_tau = mod(z_tau, (2*pi));
            tau_step_count = 0; 
        end
        
        % Rescaling step to prevent folding or overlfow etc.
        err1 = [dp1, ds1, dz1, dz_tau1]; 
        unit_err1 = err1 / sqrt(err1(1)^2 + err1(2)^2 + err1(3)^2 + err1(4)^2);
        
        err2 = [dp2, ds2, dz2, dz_tau2];
        err2 = err2 - (err2(1)*unit_err1(1) + err2(2)*unit_err1(2) + err2(3)*unit_err1(3) +...
               err2(4)*unit_err1(4))*unit_err1;
        unit_err2 = err2 / sqrt(err2(1)^2 + err2(2)^2 + err2(3)^2 + err2(4)^2);
        
        err3 = [dp3, ds3, dz3, dz_tau3];
        err3 = err3 - (err3(1)*unit_err2(1) + err3(2)*unit_err2(2) +...
                 err3(3)*unit_err2(3) + err3(4)*unit_err2(4))*unit_err2 -...
                 (err3(1)*unit_err1(1) + err3(2)*unit_err1(2) + err3(3)*unit_err1(3) +...
                 err3(4)*unit_err1(4))*unit_err1;
        unit_err3 = err3 / sqrt(err3(1)^2 + err3(2)^2 + err3(3)^2 + err3(4)^2);
        
        err4 = [dp4, dv4, dz4, dz_tau4];
        err4 = err4 - (err4(1)*unit_err3(1) + err4(2)*unit_err3(2) + err4(3)*unit_err3(3) +...
               err4(4)*unit_err3(4))*unit_err3 - (err4(1)*unit_err2(1) + err4(2)*unit_err2(2) +...
               err4(3)*unit_err2(3) + err4(4)*unit_err2(4))*unit_err2 - (err4(1)*unit_err1(1) +...
               err4(2)*unit_err1(2) + err4(3)*unit_err1(3) + err4(4)*unit_err1(4))*unit_err1;
        unit_err4 = err4 / sqrt(err4(1)^2 + err4(2)^2 + err4(3)^2 + err4(4)^2);
        
        new_vals1 = unit_err1 * err_size1;
        dp1 = new_vals1(1); ds1 = new_vals1(2); dz1 = new_vals1(3); dz_tau1 = new_vals1(4);
        
        new_vals2 = unit_err2 * err_size2;
        dp2 = new_vals2(1); ds2 = new_vals2(2); dz2 = new_vals2(3); dz_tau2 = new_vals2(4);
        
        new_vals3 = unit_err3 * err_siez3;
        dp3 = new_vals3(1); ds3 = new_vals3(2); dz3 = new_vals3(3); dz_tau3 = new_vals3(4);
        
        new_vals4 = unit_err4 * err_siez4;
        dp4 = new_vals4(1); dv4 = new_vals4(2); dz4 = new_vals3(3); dz_tau4 = new_vals4(4);
        
    end
    
    % Initialising sums for lyapunov exponents
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;
    
    % Iterations to actually calculate lyapunov exponent
    for i = 1 : num_iter
        
        % Each of the 4 steps of RK4 on each of 5 systems
        %---------------------------------------------------------------
        k1 = s;
        l1 = a*s + b*p + c*p^3 + force1*cos(z) + force2*cos(z_tau);
  
        m1 = ds1;
        n1 = a*ds1 + b*dp1 + 3*c*p^2*dp1 - force1*sin(z)*dz1 - force2*sin(z_tau)*dz_tau1;
        
        f1 = ds2;
        g1 = a*ds2 + b*dp2 + 3*c*p^2*dp2 - force1*sin(z)*dz2 - force2*sin(z_tau)*dz_tau2;
        
        h1 = ds3;
        j1 = a*ds3 + b*dp3 + 3*c*p^2*dp3 - force1*sin(z)*dz3 - force2*sin(z_tau)*dz_tau3;
        
        u1 = dv4;
        v1 = a*dv4 + b*dp4 + 3*c*p^2*dp4 - force1*sin(z)*dz4 - force2*sin(z_tau)*dz_tau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        k2 = s + 0.5*t_step*l1;
        l2 = a*(s + 0.5*t_step*l1) + b*(p + 0.5*t_step*k1) + c*(p + 0.5*t_step*k1)^3 +...
             force1*cos(z + 0.5*t_step) + force2*cos(z_tau + 0.5*tau_step);
        
        m2 = ds1 + 0.5*t_step*n1;
        n2 = a*(ds1 + 0.5*t_step*n1) + b*(dp1 + 0.5*t_step*m1) +...
             3*c*(p + 0.5*t_step*k1)^2 * (dp1 + 0.5*t_step*m1) -...
             force1*sin(z + 0.5*t_step)*dz1 - force2*sin(z_tau + 0.5*tau_step)*dz_tau1;
        
        f2 = ds2 + 0.5*t_step*g1;
        g2 = a*(ds2 + 0.5*t_step*g1) + b*(dp2 + 0.5*t_step*f1) +...
             3*c*(p + (0.5*t_step*k1))^2 * (dp2 + 0.5*t_step*f1) -...
             force1*sin(z + 0.5*t_step)*dz2 - force2*sin(z_tau + 0.5*tau_step)*dz_tau2;
        
        h2 = ds3 + (0.5*t_step*j1);
        j2 = a*(ds3 + 0.5*t_step*j1) + b*(dp3 + 0.5*t_step*h1) +...
             3*c*(p + 0.5*t_step*k1)*(p + 0.5*t_step*k1)*(dp3 + 0.5*t_step*h1) -...
             force1*sin(z + 0.5*t_step)*dz3 - force2*sin(z_tau + 0.5*tau_step)*dz_tau3;
        
        u2 = dv4 + 0.5*t_step*v1;
        v2 = a*(dv4 + 0.5*t_step*v1) + b*(dp4 + 0.5*t_step*u1) +...
             3*c*(p + 0.5*t_step*k1)^2 * (dp4 + 0.5*t_step*u1) -...
             force1*sin(z + 0.5*t_step)*dz4 - force2*sin(z_tau + 0.5*tau_step)*dz_tau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        k3 = s + 0.5*t_step*l2;
        l3 = a*(s + 0.5*t_step*l2) + b*(p + 0.5*t_step*k2) + c*(p + 0.5*t_step*k2)^3 +...
             force1*cos(z + 0.5*t_step) + force2*cos(z_tau + 0.5*tau_step);
        
        m3 = ds1 + 0.5*t_step*n2;
        n3 = a*(ds1 + 0.5*t_step*n2) + b*(dp1 + 0.5*t_step*m2) +...
             3*c*(p + 0.5*t_step*k2)^2 * (dp1 + 0.5*t_step*m2) -...
             force1*sin(z + 0.5*t_step)*dz1 - force2*sin(z_tau + 0.5*tau_step)*dz_tau1;
        
        f3 = ds2 + 0.5*t_step*g2;
        g3 = a*(ds2 + 0.5*t_step*g2) + b*(dp2 + 0.5*t_step*f2) +...
             3*c*(p + 0.5*t_step*k2)*(p + 0.5*t_step*k2)*(dp2 + 0.5*t_step*f2) -...
             force1*sin(z + 0.5*t_step)*dz2 - force2*sin(z_tau + 0.5*tau_step)*dz_tau2;
        
        h3 = ds3 + 0.5*t_step*j2;
        j3 = a*(ds3 + 0.5*t_step*j2) + b*(dp3 + 0.5*t_step*h2) +...
            3*c*(p + 0.5*t_step*k2)*(p + 0.5*t_step*k2)*(dp3 + 0.5*t_step*h2) -...
            force1*sin(z + 0.5*t_step)*dz3 - force2*sin(z_tau + 0.5*tau_step)*dz_tau3;
        
        u3 = dv4 + 0.5*t_step*v2;
        v3 = a*(dv4 + 0.5*t_step*v2) + b*(dp4 + 0.5*t_step*u2) +...
             3*c*(p + 0.5*t_step*k2)*(p + 0.5*t_step*k2)*(dp4 + 0.5*t_step*u2) -...
             force1*sin(z + 0.5*t_step)*dz4 - force2*sin(z_tau + 0.5*tau_step)*dz_tau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        k4 = s + t_step*l3;       
        l4 = a*(s + t_step*l3) + b*(p + t_step*k3) + c*(p + t_step*k3)^3 +...
             force1*cos(z + t_step) + force2*cos(z_tau + tau_step);
       
        m4 = ds1 + t_step*n3;
        n4 = a*(ds1 + t_step*n3) + b*(dp1 + t_step*m3) +...
            3*c*(p + t_step*k3)^2 * (dp1 + t_step*m3) -...
            force1*sin(z + t_step)*dz1 - force2*sin(z_tau + tau_step)*dz_tau1;
        
        f4 = ds2 + t_step*g3;
        g4 = a*(ds2 + t_step*g3) + b*(dp2 + t_step*f3) +...
             3*c*(p + t_step*k3)^2 * (dp2 + (t_step*f3)) -...
             force1*sin(z + t_step)*dz2 - force2*sin(z_tau + tau_step)*dz_tau2;
        
        h4 = ds3 + t_step*j3;
        j4 = a*(ds3 + t_step*j3) + b*(dp3 + t_step*h3) +...
            3*c*(p + t_step*k3)^2 * (dp3 + t_step*h3) -...
            force1*sin(z + t_step)*dz3 - force2*sin(z_tau + tau_step)*dz_tau3;
        
        u4 = dv4 + t_step*v3;
        v4 = a*(dv4 + t_step*v3) + b*(dp4 + t_step*u3) +...
             3*c*(p + t_step*k3)*(p + t_step*k3)*(dp4 + t_step*u3) -...
             force1*sin(z + t_step)*dz4 - force2*sin(z_tau + tau_step)*dz_tau4;
        %---------------------------------------------------------------
        
        p = p + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*t_step;
        s = s + (1/6)*(l1 + 2*l2 + 2*l3 + l4)*t_step;
        
        dp1 = dp1 + (1/6)*(m1 + 2*m2 + 2*m3 + m4)*t_step;
        ds1 = ds1 + (1/6)*(n1 + 2*n2 + 2*n3 + n4)*t_step;

        dp2 = dp2 + (1/6)*(force1 + 2*f2 + 2*f3 + f4)*t_step;
        ds2 = ds2 + (1/6)*(g1 + 2*g2 + 2*g3 + g4)*t_step;

        dp3 = dp3 + (1/6)*(h1 + 2*h2 + 2*h3 + h4)*t_step;
        ds3 = ds3 + (1/6)*(J1 + 2*j2 + 2*j3 + j4)*t_step;
        
        dp4 = dp4 + (1/6)*(u1 + 2*u2 + 2*u3 + u4)*t_step;
        dv4 = dv4 + (1/6)*(v1 + 2*v2 + 2*v3 + v4)*t_step;
        
        z = z + t_step;
        t_step_count = t_step_count + 1;
        
        z_tau = z_tau + tau_step;
        tau_step_count = tau_step_count + 1;
        
        if t_step_count == step_size
            z = 0;
            t_step_count = 0;
        end
        
        if tau_step_count == step_size
            z_tau = mod(z_tau, (2*pi));
            tau_step_count = 0; 
        end
        
        % Rescaling step to prevent folding or overlfow etc.
        err1 = [dp1, ds1, dz1, dz_tau1]; 
        unit_err1 = err1 / sqrt(err1(1)^2 + err1(2)^2 + err1(3)^2 + err1(4)^2);
        
        err2 = [dp2, ds2, dz2, dz_tau2];
        err2 = err2 - (err2(1)*unit_err1(1) + err2(2)*unit_err1(2) + err2(3)*unit_err1(3) +...
               err2(4)*unit_err1(4))*unit_err1;
        unit_err2 = err2 / sqrt(err2(1)^2 + err2(2)^2 + err2(3)^2 + err2(4)^2);
        
        err3 = [dp3, ds3, dz3, dz_tau3];
        err3 = err3 - (err3(1)*unit_err2(1) + err3(2)*unit_err2(2) +...
                 err3(3)*unit_err2(3) + err3(4)*unit_err2(4))*unit_err2 -...
                 (err3(1)*unit_err1(1) + err3(2)*unit_err1(2) + err3(3)*unit_err1(3) +...
                 err3(4)*unit_err1(4))*unit_err1;
        unit_err3 = err3 / sqrt(err3(1)^2 + err3(2)^2 + err3(3)^2 + err3(4)^2);
        
        err4 = [dp4, dv4, dz4, dz_tau4];
        err4 = err4 - (err4(1)*unit_err3(1) + err4(2)*unit_err3(2) + err4(3)*unit_err3(3) +...
               err4(4)*unit_err3(4))*unit_err3 - (err4(1)*unit_err2(1) + err4(2)*unit_err2(2) +...
               err4(3)*unit_err2(3) + err4(4)*unit_err2(4))*unit_err2 - (err4(1)*unit_err1(1) +...
               err4(2)*unit_err1(2) + err4(3)*unit_err1(3) + err4(4)*unit_err1(4))*unit_err1;
        unit_err4 = err4 / sqrt(err4(1)^2 + err4(2)^2 + err4(3)^2 + err4(4)^2);
        
        new_vals1 = unit_err1 * err_size1;
        dp1 = new_vals1(1); ds1 = new_vals1(2); dz1 = new_vals1(3); dz_tau1 = new_vals1(4);
        
        new_vals2 = unit_err2 * err_size2;
        dp2 = new_vals2(1); ds2 = new_vals2(2); dz2 = new_vals2(3); dz_tau2 = new_vals2(4);
        
        new_vals3 = unit_err3 * err_siez3;
        dp3 = new_vals3(1); ds3 = new_vals3(2); dz3 = new_vals3(3); dz_tau3 = new_vals3(4);
        
        new_vals4 = unit_err4 * err_siez4;
        dp4 = new_vals4(1); dv4 = new_vals4(2); dz4 = new_vals3(3); dz_tau4 = new_vals4(4);
        
        % Adding to sum to calculate lyapunov exponents.
        sum1 = sum1 + log((sqrt(err1(1)^2 + err1(2)^2 + err1(3)^2 + err1(4)^2)) / err_size1);
        sum2 = sum2 + log((sqrt(err2(1)^2 + err2(2)^2 + err2(3)^2 + err2(4)^2)) / err_size2);
        sum3 = sum3 + log((sqrt(err3(1)^2 + err3(2)^2 + err3(3)^2 + err3(4)^2)) / err_siez3);
        sum4 = sum4 + log((sqrt(err4(1)^2 + err4(2)^2 + err4(3)^2 + err4(4)^2)) / err_siez4);
        
    end 
    
    % Taking the average to get final answer for LE and put into array
    lyap_exponent1 = sum1 / (num_iter*t_step);
    lyap_exponent2 = sum2 / (num_iter*t_step);
    lyap_exponent3 = sum3 / (num_iter*t_step);
    lyap_exponent4 = sum4 / (num_iter*t_step);
    
    lyap_exponents = [lyap_exponent1, lyap_exponent2, lyap_exponent3, lyap_exponent4];
    
end