% Function to plot the poincare section of the oscillator, or the
% stroboscopic view as you may call it. Performs numerical integration and
% records the solution after each cycle of 2*pi, but can start at time != 0
% to take a different "slice" of the oscillators orbit, chosen by the input
% value of phase irrational frequency is an input here, usually golden
% ratio but e and sqrt(7) / 7 produce nice results too. steps per cycles set to
% 50 for near optimum speed/accuracy ratio.

% !ATTENTION! - For wolfram team: some code here requires a background
% knowledge of the quasi-periodically driven duffing oscillator. Here, the
% parameters of the ODE system you are finding the poincare setion for are
% described below, and the irrational frequenct or irrationalFreq is the
% drving freqiency to generate quasi-periodic motion. the input value
% "numCycles" is how many driving cycles of t = 2*pi you want for the forcing. 

% d2y/dx2 = a(dy/dx) + bx + cx^3 + d1Cos(t) + d2Cos(tau)
% suggested: a = -0.2, b = 1, c = -1, d1 = 0.3, d2 = 0.3, irrationalFreq =
% (sqrt(5) - 1) / 2, numCycles = 10^5


function [position, speed, tauVector] = PoincareByParameters (a, b, c, d1, d2, irrationalFreq, numCycles)
    
    t = 0;
    tau = t * irrationalFreq;
    stepsPerCycle = 50;
  
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
    plot3(tauVector((0.02 * numCycles) : numCycles), position((0.02 * numCycles) : numCycles),...
        speed((0.02 * numCycles) : numCycles), '+', 'MarkerSize', 0.1);
%     plot(tauVector((0.001*numCycles) : numCycles), ...
%         position((0.001*numCycles) : numCycles), '+', 'MarkerSize', 0.1);
% 
    xlabel('tau (modified time)');
    ylabel('x (Position)');
    zlabel('dx/dt (speed)');
    title(['Poincare section at d2 = ' num2str(d2) '.']);
    
end

% Function to calculate the correlation dimension of the attractor. Inputs
% are the parameters of system, and "type" which is what type of attractor
% it should be. 1 = x vs tau, 2 = x vs dx/dt, 3 = all three (3D plot).
% output is the fractal dimension itself.

% Function works by measuring "density" or probability that points will be
% near each other, within a given radius. This is done by setting a radius,
% looking at each point on the attractor and seeing how many other points
% are within a radius from it. This number is continually summed as we
% observe each point in the attractor, and then divided by number of points
% * number of points - 1 to give an approximation to the correlation
% integral. The dimension is then found from the logs of these values
% against the logs of the size of radii. radius size goes from 0.005 to
% 0.0025, any smaller and the attractor simply isnt precise enough for
% calculations.


function [fractalDim] = CorrelationDimension (a, b, c, d1, d2, irrationalFreq, type)
    
    % number of iterations with varying radii
    numTests = 25;
    radius   = linspace(0.005, 0.0025, numTests);
    
    [position, speed, tauVector] = PoincareByParameters(a, b, c, d1, d2, irrationalFreq, 10^4);
    numPoints = length(position);
    
    corrSums = zeros(1, numTests);
    
    position  = position((0.02 * numPoints) : numPoints);
    speed     = speed((0.02 * numPoints) : numPoints);
    tauVector = tauVector((0.02 * numPoints) : numPoints);
    
    numPoints = length(position);
    
    % 3 = all three dimensions included
    if type == 3
        
        for i = 1 : numTests

            r = radius(i);
            numInRadius = 0;
            
            % cheking every other sample of attractor
            for j = 1 : numPoints

                p1 = position(j);
                s1 = speed(j);
                t1 = tauVector(j);
                
                % setting up the rest of the points as a vector...
                p2 = position(j + 1 : numPoints);
                s2 = speed(j + 1 : numPoints);
                t2 = tauVector(j + 1 : numPoints);
                
                % ...So it saves a lot of time by running in parallel.
                distance = ((p1 - p2) .* (p1 - p2)) + ...
                ((s1 - s2) .* (s1 - s2)) + ((t1 - t2) .* (t1 - t2));

                % array of bits that mimicks heaviside step function
                isInRadius  = (r^2 - distance) >= 0;
                numInRadius = numInRadius + sum(isInRadius);

            end 
            
            % correlation sum, approximation to correlation intgral.
            corrSums(i) = log((2 * numInRadius) / (numPoints * (numPoints - 1)));
            
            plot3(position, speed, tauVector, '.', 'MarkerSize', 1);
            xlabel('x');
            ylabel('dx/dt');
            zlabel('tau');
            
        end
        
    % 2 = x vs dx/dt, see above for more details.
    elseif type == 2
        
        for i = 1 : numTests

            r = radius(i);
            numInRadius = 0;

            for j = 1 : numPoints

                p1   = position(j);
                s1   = speed(j);

                p2   = position(j + 1 : numPoints);
                s2   = speed(j + 1 : numPoints);

                distance = ((p1 - p2) .* (p1 - p2)) + ...
                ((s1 - s2) .* (s1 - s2));

                isInRadius = (r^2 - distance) > 0;
                numInRadius = numInRadius + sum(isInRadius);

            end 
            
            corrSums(i) = log((2 * numInRadius) / (numPoints * (numPoints - 1)));
            
            plot(position, speed, '.', 'MarkerSize', 1);
            xlabel('x');
            ylabel('dx/dt');

        end
        
    % 1 = x vs tau.
    elseif type == 1
        
        for i = 1 : numTests

            r = radius(i);
            numInRadius = 0;

            for j = 1 : numPoints

                p1 = position(j);
                t1 = tauVector(j);

                p2 = position(j + 1 : numPoints);
                t2 = tauVector(j + 1 : numPoints);

                distance = ((p1 - p2) .* (p1 - p2)) + ...
                    ((t1 - t2) .* (t1 - t2));

                isInRadius  = (r^2 - distance) > 0;
                numInRadius = numInRadius + sum(isInRadius);

            end 

            corrSums(i) = log((2 * numInRadius) / (numPoints * (numPoints - 1)));
            
            plot(tauVector, position, '.', 'MarkerSize', 0.1);
            xlabel('tau');
            ylabel('x');

        end
        
    % only 1, 2, 3 are accepted attractor space dimesnions.
    else
        disp('enter correct dimension of attractor space')
        return
    end
 
    % setting up recorded values for linear regression
    corrSums = transpose(corrSums);
    radius   = transpose(radius);
    
    % fractal dimension as gradient via basic linear regression
    fractalDim = [ones(numTests, 1), log(radius)] \ corrSums;
    fractalDim = fractalDim(2);
    
    title(['Correlation dimension = ' num2str(fractalDim) '.']);
     
end