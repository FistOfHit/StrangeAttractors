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
% against the logs of the size of radii.

function [fractalDim] = CorrelationDimension (a, b, c, d1, d2, irrationalFreq, numCycles)
    
    % number of iterations for varying radii
    numTests = 25;
    
    % Generating data and removing transient effects
    [position, speed, tauVector] = PoincareByParameters(a, b, c, d1, d2, irrationalFreq, numCycles);
    
    position  = position((0.2 * numCycles) : numCycles);
    speed     = speed((0.2 * numCycles) : numCycles);
    tauVector = tauVector((0.2 * numCycles) : numCycles);
    numPoints = length(position);
    
    corrSums = zeros(1, numTests);
    radii   = 10 .^ linspace(-1.5, -2.5, numTests);

    for i = 1 : numTests

        radius = radii(i);
        numInRadius = 0;

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
            isInRadius  = ((radius ^ 2) - distance) > 0 ;
            numInRadius = numInRadius + sum(isInRadius);

        end 

        % correlation sum, approximation to correlation intgral.
        corrSums(i) = (1 * numInRadius) / (numPoints * (numPoints  - 1));

%         plot3(position, speed, tauVector, '.', 'MarkerSize', 1);
%         xlabel('x');
%         ylabel('dx/dt');
%         zlabel('tau');

    end
   
    % setting up recorded values for linear regression
    corrSums = log(corrSums);
    radii    = log(radii);
    
    % fractal dimension as gradient via basic linear regression
    fractalDim = polyfit(radii, corrSums, 1);
    fractalDim = fractalDim(1);
    
%     title(['Correlation dimension = ' num2str(fractalDim) '.']);
    
%     plot(radii, corrSums);
%     hold on
%     plot(radii, (fractalDim(1) .* radii) + fractalDim(2), '+');
%     hold off

end