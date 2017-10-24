% Function to calculate the "sphere count dimension", using spheres to
% cover the points on the attractor rather than boxes. Function takes
% parameters of the oscillator, as well as another paramter called "type"
% which determines the type of attractor you get. 1 = x vs tau, 2 = x vs
% dx/dt, 3 = all three of them (3D plot). Output is the sphere counting
% dimension (Hausdorff dimension?) of the attractor, also included in the
% title of each plot.

% The use of spheres rather than boxes simplifies calculations and reduces
% error, win-win. Algorithm works by putting a sphere around a point,
% listing the points inside the sphere as "covered" and then doing the same
% for other points that are not covered as such. This coveres the attractor
% with spheres which do overlap. The dimension is then calculated from the
% number of spheres and their radii. more details and technicalities will
% be explained below.

function [fractalDim] = SphereCountDimension(a, b, c, d1, d2, irrationalFreq, numCycles)

    [position, speed, tauVector] = PoincareByParameters (a, b, c, d1, d2, irrationalFreq, numCycles);

    numPoints = length(position);
    % removing first 20% of points in attractor, removing transient effects
    position  = position((0.2 * numPoints) : numPoints);
    speed     = speed((0.2 * numPoints) : numPoints);
    tauVector = tauVector((0.2 * numPoints) : numPoints);
    numPoints = length(position);
    
    numTests = 25;

    radii = 10 .^ linspace(0, -1.3, numTests);
    spheresNeeded  = zeros(1, numTests);
    isPointCovered = zeros(1, numPoints);
    numSpheres = 0;

    for i = 1 : numTests

        % going through each point in attractor
        for j = 1 : numPoints

            p_i = position(j);
            s_i = speed(j);
            t_i = tauVector(j);

            radius = radii(i);

            if isPointCovered(j) == 0

                % if not covered, then cover it
                isPointCovered(j) = 1;
                numSpheres = numSpheres + 1;

                % check what other points are in this sphere we made
                p_j = position(j + 1 : numPoints);
                s_j = speed(j + 1 : numPoints);
                t_j = tauVector(j + 1 : numPoints);

                coverage = (((p_i - p_j).^2) + ((s_i - s_j).^2)  + ((t_i - t_j).^2))  < radius^2;
                isPointCovered(j + 1 : numPoints) = isPointCovered(j + 1 : numPoints) + coverage;

                if sum(isPointCovered(j + 1 : numPoints) > 1) ~= 0

                    isPointCovered(j + 1 : numPoints) = isPointCovered(j + 1 : numPoints) - coverage;

                    numSpheres = numSpheres - 1;
                    isPointCovered(j) = 0;

                end

            end

        end

        spheresNeeded(i) = numSpheres;
        numSpheres = 0;

        isPointCovered = zeros(1, numPoints);

    end

%     % Here is your attractor for visiual inspection!
%     plot3(position, speed, tauVector, '.', 'MarkerSize', 1);
%     xlabel('x');
%     ylabel('dx/dt');
%     zlabel('tau');

    % setting up the values found for linear regression
    xAxis = log(1 ./ radii);
    yAxis = log(spheresNeeded);

    % gradient (fractal dimension) obtained by linear regression
    fractalDim = polyfit(xAxis, yAxis, 1);
    fractalDim = fractalDim(1);
    
%     plot(xAxis, yAxis);
%     hold on
%     plot(xAxis, fractalDim .* xAxis, '+');
%     hold off

end



