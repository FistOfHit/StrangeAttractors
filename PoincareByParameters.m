% Function to plot the poincare section of the oscillator, or the
% stroboscopic view as you may call it. Performs numerical integration
% and records the solution after each cycle of 2*pi. Outputs data in 3 
% 1xnumCycles sized vectors, but also can plot the data. Inputs are 
% the parameters of the oscillator, and the number of driving cycles
% needed. steps per cycles set to 25.

function [position, speed, tauVector] = PoincareByParameters (a, b, c, d1, d2, irrationalFreq, numCycles)
    
    t = 0;
    tau = t * irrationalFreq;
    stepsPerCycle = 25;
  
    h = (2*pi) / stepsPerCycle;
    hTau = h * irrationalFreq;
    
    numIterations = stepsPerCycle * numCycles;
    
    position  = zeros(1, numCycles);
    speed     = zeros(1, numCycles);
    tauVector = zeros(1, numCycles);
    
    % counter to kep track of progression through driving cycle
    hCount = 0;
    hTauCount = 0;
    
    % Initial conditons  
    p = 1;
    s = 0;
    
    for i = 1 : numIterations

        % Runge-Kutta calculations (Hard-coded for performance)----------
        K1 = s;
        L1 = a*s + b*p + c*p*p*p + d1*cos(t) + d2*cos(tau);

        K2 = s + (0.5*h*L1);
        L2 = a*(s + (0.5*h*L1)) + b*(p + (0.5*h*K1)) +...
            c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(p + (0.5*h*K1)) +...
            d1*cos(t + (0.5*h)) + d2*cos(tau + (0.5*hTau));

        K3 = s + (0.5*h*L2);
        L3 = a*(s + (0.5*h*L2)) + b*(p + (0.5*h*K2)) +...
            c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(p + (0.5*h*K2)) +...
            d1*cos(t + (0.5*h)) + d2*cos(tau + (0.5*hTau));

        K4 = s + (h*L3);       
        L4 = a*(s + (h*L3)) + b*(p + (h*K3)) +...
            c*(p + (h*K3))*(p + (h*K3))*(p + (h*K3)) +...
            d1*cos(t + h) + d2*cos(tau + hTau);

        p = p + (1/6)*(K1 + (2*K2) + (2*K3) + K4)*h;
        s = s + (1/6)*(L1 + (2*L2) + (2*L3) + L4)*h;
        % ----------------------------------------------------------------
        
        t = t + h;
        hCount = hCount + 1;
        
        tau = tau + hTau;
        hTauCount = hTauCount + 1;
        
        if hCount == stepsPerCycle
            
            t = 0;
            hCount = 0;
            
            % recording solution after every cycle
            position((i / stepsPerCycle)) = p;
            speed((i / stepsPerCycle)) = s;
            
        end
        
        if hTauCount == stepsPerCycle 
            
           tau = mod(tau, (2*pi));
           hTauCount = 0; 
           
           % recording solution after every cycle
           tauVector((i / stepsPerCycle)) = tau;
           
        end

    end
    
%     % Plot in either 3D or 2D depending on need
    plot3(position((0.2 * numCycles) : numCycles), speed((0.2 * numCycles) : numCycles),...
        tauVector((0.2 * numCycles) : numCycles), '.', 'MarkerSize', 0.1);
%     plot(position((0.2*numCycles) : numCycles), ...
%         speed((0.2*numCycles) : numCycles), '+', 'MarkerSize', 0.1);

    xlabel('x');
    ylabel('dx/dt');
    zlabel('tau');
    title('Poincare section at t = 2pi');
    
end