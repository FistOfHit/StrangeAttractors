% Function to solve the Quasi-periodically driven duffing oscillator with
% 4th order Runge-Kutta numerical integration. Outputs a vector of the x
% values and dx/dt values, and takes in the number of steps per a cycle
% (every 2*pi time units) which is effectively the precision, number of
% 2*pi time units you want to solve for and the initial conditions as a 1x2
% vector. Utilises two time dimensions, t and tau, resetting both when they
% reach 2*pi and 2*pi*GoldenRatio respectively.

function[position, speed] = DuffingSolution(stepsPerCycle, numCycles, initConditions)
    
    t   = 0;
    tau = 0;
    
    numIterations = stepsPerCycle * numCycles;
    
    position = zeros(1, (numIterations + 1));
    speed    = zeros(1, (numIterations + 1));
    
    position(1) = initConditions(1);
    speed(1)    = initConditions(2);
                        
    % Parameters of the oscillator (dy/dt = ay + bx + cx^3 + d1cos(t) +
    % d2cos(tau)) and the irrational driving frequency.
    a = -0.2; b = 1; c = -2; d1 = 0.3; d2 = 0; irrationalFreq = ((sqrt(5) - 1) / 2);

    % time steps for both time dimensions, with counters to keep track
    tStep = (2 * pi) / stepsPerCycle;
    tauStep = tStep * irrationalFreq;
    
    hCounter = 0;
    hTauCounter = 0;
    
    for i = 1 : numIterations
        
        p = position(i);
        s = speed(i);
    
        % ... please understand.
        k1 = s;
        l1 = a*s + b*p + c*p*p*p + d1*cos(t) + d2*cos(tau);

        k2 = s + (0.5*tStep*l1);
        l2 = a*(s + (0.5*tStep*l1)) + b*(p + (0.5*tStep*k1)) +...
            c*(p + (0.5*tStep*k1))*(p + (0.5*tStep*k1))*(p + (0.5*tStep*k1)) +...
            d1*cos(t + (0.5*tStep)) + d2*cos(tau + (0.5*tauStep));

        k3 = s + (0.5*tStep*l2);
        l3 = a*(s + (0.5*tStep*l2)) + b*(p + (0.5*tStep*k2)) +...
            c*(p + (0.5*tStep*k2))*(p + (0.5*tStep*k2))*(p + (0.5*tStep*k2)) +...
            d1*cos(t + (0.5*tStep)) + d2*cos(tau + (0.5*tauStep));

        k4 = s + (tStep*l3);       
        l4 = a*(s + (tStep*l3)) + b*(p + (tStep*k3)) +...
            c*(p + (tStep*k3))*(p + (tStep*k3))*(p + (tStep*k3)) +...
            d1*cos(t + tStep) + d2*cos(tau + tauStep);

        position(i + 1) = p + (1/6)*(k1 + (2*k2) + (2*k3) + k4)*tStep;
        speed(i + 1)    = s + (1/6)*(l1 + (2*l2) + (2*l3) + l4)*tStep;

        t = t + tStep;
        hCounter = hCounter + 1;
        
        tau = tau + tauStep;
        hTauCounter = hTauCounter + 1;
        
        if hCounter == stepsPerCycle
            t = 0;
            hCounter = 0;
        end
        
        if hTauCounter == stepsPerCycle 
            tau = mod(tau, (2 * pi));
            hTauCounter = 0;   
        end
        
    end
    
end





