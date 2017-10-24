% Function to plot the poincare section of the oscillator, or the
% stroboscopic view as you may call it. Performs numerical integration and
% records the solution after each cycle of 2*pi, but can start at time != 0
% to take a different "slice" of the oscillators orbit, chosen by the input
% value of phase irrational frequency is an input here, usually golden
% ratio but e and sqrt(7) produce nice results too. steps per cycles set to
% 25 for near optimum speed/accuracy ratio.

function[] = PoincareSectionPlot(numCycles, irrationalFreq, phase)

    t = phase;
    tau = t * irrationalFreq;
    stepSize = 50;
  
    h = (2*pi) / stepSize;
    hTau = h * irrationalFreq;
    
    numIterations = stepSize * numCycles;
    
    position  = zeros(1, numCycles);
    speed     = zeros(1, numCycles);
    tauVector = zeros(1, numCycles);
         
    % Parameters of the oscillator (dy/dt = ay + bx + cx^3 + d1cos(t) +
    % d2cos(tau))
    a = -0.05; b = 0; c = -1; d1 = 7.5; d2 = 0;
    
    % counter to kep track of progression through driving cycle
    hCount = 0;
    hTauCount = 0;
    
    % Initial conditons  
    p = 1;
    s = 0;
    
    for i = 1 : numIterations

        % Runge-Kutta calculations---------------------------------------
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
        
        if hCount == stepSize
            
            t = phase;
            hCount = 0;
            
            % recording solution after every cycle
            position((i / stepSize)) = p;
            speed((i / stepSize)) = s;
            
        end
        
        if hTauCount == stepSize 
            
           tau = mod(tau, (2*pi));
           hTauCount = 0; 
           
           % recording solution after every cycle
           tauVector((i / stepSize)) = tau;
           
        end

    end
    
    % Plot in either 3D or 2D depending on need
    plot3(tauVector((0.2*numCycles) : numCycles), position((0.2*numCycles) : numCycles),...
        speed((0.2*numCycles) : numCycles), '+', 'MarkerSize', 0.00000000001);
%     plot(tauVector((0.1*numCycles) : numCycles), ...
%         position((0.1*numCycles) : numCycles), '+', 'MarkerSize', 0.00000000001);

%         plot(position((0.1*numCycles) : numCycles), ...
%         speed((0.1*numCycles) : numCycles), '+', 'MarkerSize', 0.00000000001);

    xlabel('tau (modified time)');
    ylabel('x (Position)');
    zlabel('dx/dt (speed)');
    title(['Poincare section at t = ' num2str(phase / pi) 'pi.']);
    
end