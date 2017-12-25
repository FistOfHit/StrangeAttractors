function[pos, speed, tau_array] = Create_attractor(a, b, c, force1, force2, irr_freq, num_cycles)
% Generate attractor set of duffing oscillator from t = 0 to 2*pi*num_cycles.
%
% Inputs:
% a = -1*Damping strength
% b = Linear stiffness
% c = non-linearity in restoring force
% force1 = amplitude of force in rational frequency driving force (2*pi Hz)
% force2 = amplitude of force in irrational frequecy driving force (2*pi*irr_freq) Hz)
% irr_freq = value of irrational driving frequency 
% num cycles = number of 2*pi long driving cycles of the forcing function
%
% Outputs:
% pos = position - array of x-values after each cycle
% speed = array of dx/dt values after each cycle
% tau_array = array of tau values after each cycle
%
% Figures:
% Plot of attractor set points in 3D
%
% Here the values are generated from the Duffing_solution function with the parameters given, and
% are plotted in either 2 or 3D depending on need. For more information on how the solution is
% found, check out either the documentation in the repository or the description of the
% Duffing_solution function. Here, we only plot the solution after each reset of the two temporal
% dimensions, to simulate taking a Poincare section. The final set is plotted after remving the
% transients ro show long term behaviour, giving an attractor set.

    steps_per_cyc = 100;
    t_step = (2*pi)/steps_per_cyc;
    tau_step = t_step*irr_freq;
    
    num_iter = steps_per_cyc * num_cycles;
    
    pos = zeros(1, num_cycles);
    speed = zeros(1, num_cycles);
    tau_array = zeros(1, num_cycles);
    
    % counter to kep track of progression through driving cycle
    t_step_count = 0;
    tau_step_count = 0;
    
    % Initial conditons (default at x = 1, dx/dt = 0)
    p = 1;
    s = 0;
    t = 0;
    tau = 0;
    
    for i = 1 : num_iter

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
        
        t = t + t_step;
        t_step_count = t_step_count + 1;
        
        tau = tau + tau_step;
        tau_step_count = tau_step_count + 1;
        
        % Recording values only after each cycle
        if t_step_count == steps_per_cyc
            t = 0;
            t_step_count = 0;
            
            pos(i/steps_per_cyc) = p;
            speed(i/steps_per_cyc) = s;
        end
        
        if tau_step_count == steps_per_cyc 
           tau = mod(tau, (2*pi));
           tau_step_count = 0; 
           
           tau_array((i / steps_per_cyc)) = tau;
        end

    end
    
    % Plotting attractor set (means transients must be removed)
    plot3(pos(0.2*num_cycles : num_cycles), speed(0.2*num_cycles : num_cycles),...
              tau_array(0.2*num_cycles : num_cycles), '.', 'MarkerSize', 0.1);
    xlabel('x');
    ylabel('dx/dt');
    zlabel('tau');
    title('Poincare section at t = 2pi');
    
end