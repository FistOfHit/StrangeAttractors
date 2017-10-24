% Function allowing you to find solution to system with the parameters as
% inputs, allowing for much greater convnience when you want to observe the
% effects on the solutions by variation of specfic or all parameters. steps
% per cycle is set to 100, and initial conditons are 1 and 0.

function [position, speed] = DuffingSolutionByParameters (a, b, c, d1, d2, irrationalFreq, initConditions)

    t = 0;
    tau = 0;
    
    stepsPerCycle = 100;
    numCycles = (10^4);
    h = (2*pi) / stepsPerCycle;
    hTau = h * irrationalFreq;
    
    numIterations = stepsPerCycle * numCycles;
   
    position = zeros(1, (numIterations + 1));
    speed    = zeros(1, (numIterations + 1));
    
    % Initial conditons  
    position(1) = initConditions(1);
    speed(1)    = initConditions(2);
                        
    % Counters to track progresion through driving cycles
    hCount = 0;
    hTauCount = 0;
    
    for i = 1 : numIterations
        
        p = position(i);
        s = speed(i);

        k1 = s;
        l1 = a*s + b*p + c*p*p*p + d1*cos(t) + d2*cos(tau);

        k2 = s + (0.5*h*l1);
        l2 = a*(s + (0.5*h*l1)) + b*(p + (0.5*h*k1)) +...
            c*(p + (0.5*h*k1))*(p + (0.5*h*k1))*(p + (0.5*h*k1)) +...
            d1*cos(t + (0.5*h)) + d2*cos(tau + (0.5*hTau));

        k3 = s + (0.5*h*l2);
        l3 = a*(s + (0.5*h*l2)) + b*(p + (0.5*h*k2)) +...
            c*(p + (0.5*h*k2))*(p + (0.5*h*k2))*(p + (0.5*h*k2)) +...
            d1*cos(t + (0.5*h)) + d2*cos(tau + (0.5*hTau));

        k4 = s + (h*l3);       
        l4 = a*(s + (h*l3)) + b*(p + (h*k3)) +...
            c*(p + (h*k3))*(p + (h*k3))*(p + (h*k3)) +...
            d1*cos(t + h) + d2*cos(tau + hTau);

        position(i + 1) = p + (1/6)*(k1 + (2*k2) + (2*k3) + k4)*h;
        speed(i + 1) = s + (1/6)*(l1 + (2*l2) + (2*l3) + l4)*h;

        t = t + h;
        hCount = hCount + 1;
        
        tau = tau + hTau;
        hTauCount = hTauCount + 1;
        
        if hCount == stepsPerCycle
            t = 0;
            hCount = 0;
        end
        
        if hTauCount == stepsPerCycle
           tau = mod(tau, (2*pi));
           hTauCount = 0;   
        end
        
    end
    
end
