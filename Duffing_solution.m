function[pos, speed, tau_array] = Duffing_solution(a, b, c, force1, force2, irr_freq, init_cond)
% Find Solution over time of duffing oscillator with parameters and initial condtions given
%
% Inputs:
% a = -1*Damping strength
% b = Linear stiffness
% c = non-linearity in restoring force
% force1 = amplitude of force in rational frequency driving force (2*pi Hz)
% force2 = amplitude of force in irrational frequecy driving force (2*pi*irr_freq) Hz)
% irr_freq = value of irrational driving frequency 
% init_cond = initial conditions for dynamical system
% 
% Outputs:
% pos = array of solutions, x-values
% speed = array of solutions, dx/dt values
% tau_array = array of solutions, tau values
% 
% The solution is found and output for each time step in arrays corresponding to each non (standard)
% time dimension, hence position (x), speed (dx/dt) and tau dimensions. The soltuions are found from
% the parameters provided by using Runge-Kutta fourth order numerical integration on the two first 
% order ODE's split from one second order ODE. As they have inter dependancies, we must run two RK4
% methods at the same time, as denoted by the k's and j's. As discussed in greater detail in the
% full official report available in this repository, time and tau steps are counted and reset once
% the value in those dimensions reaches 2*pi with conditional statements. 
% 
% Note that here, for simplicity the Quasi-periodically driven duffing oscillator is represented as:
%                 d2x/dt2 = a*dx/dt + b*x + c*x^3 + force1*cos(t) + force2*cos(tau)
% as opposed to the more meaningful but complicated notation found in some other sources. With this,
% there are some limitations to be set on the parameters too:
% a: from 0 to -infinity only, setting a +ve will result in a non-physical system
% b: from 0 to +infintiy only, setting b -ve will result in a non-physical system
% c: from 0 to -infinity only, setting c +ve will result in a non-physical system
% force1 and force2: non-zero, else there is no driving force
% irr_freq: non-zero, else there is no quasi-periodicity occuring
% init_cond: preferably between +-10, else its just a pointless waste of time waiting for system to
% "calm down" and show any interesting results.


    % Number of time steps per each 2*pi length of time (1 driving cycle)
    steps_per_cyc = 100;
    % Number of 2*pi long cycles of the driving force 
    num_cycles = 10^4;
    num_iter = steps_per_cyc*num_cycles;
    
    pos = zeros(1, num_iter + 1);
    speed = zeros(1, num_iter + 1);
    tau_array = zeros(1, num_iter + 1);
    
    % Initial conditons  
    p = init_cond(1);
    s = init_cond(2);
    t = 0;
    tau = 0;   
    
    % Time step sizes in both time dimensions 
    t_step = (2*pi)/steps_per_cyc;
    tau_step = t_step*irr_freq;
    % Counters to track progresion through driving cycles
    t_step_count = 0;
    tau_step_count = 0;
    
    for i = 1 : num_iter
        
        pos(i) = p;
        speed(i) = s;

        % Runge Kutta 4th order numerical integration for this system (2 ODE's simultaneous)
        k1 = s;
        j1 = a*s + b*p + c*(p^3) + force1*cos(t) + force2*cos(tau);

        k2 = s + 0.5*t_step*j1;
        j2 = a*(s + 0.5*t_step*j1) + b*(p + 0.5*t_step*k1) + c*(p + 0.5*t_step*k1)^3 +...
             force1*cos(t + 0.5*t_step) + force2*cos(tau + 0.5*tau_step);

        k3 = s + 0.5*t_step*j2;
        j3 = a*(s + 0.5*t_step*j2) + b*(p + 0.5*t_step*k2) + c*(p + 0.5*t_step*k2)^3 +...
             force1*cos(t + (0.5*t_step)) + force2*cos(tau + (0.5*tau_step));

        k4 = s + t_step*j3;       
        j4 = a*(s + t_step*j3) + b*(p + t_step*k3) + c*(p + t_step*k3)^3 +...
             force1*cos(t + t_step) + force2*cos(tau + tau_step);

        p = p + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*t_step;
        s = s + (1/6)*(j1 + 2*j2 + 2*j3 + j4)*t_step;

        % Resetting time after t = 2*pi
        t = t + t_step;
        t_step_count = t_step_count + 1;
        if t_step_count == steps_per_cyc
            t = 0;
            t_step_count = 0;
        end
        
        % Resetting tau after tau = 2*pi*irr_freq
        tau = tau + tau_step;
        tau_step_count = tau_step_count + 1;
        if tau_step_count == steps_per_cyc
           tau = mod(tau, 2*pi);
           tau_array(i/steps_per_cyc) = tau;
           tau_step_count = 0;   
        end
        
    end
    
end
