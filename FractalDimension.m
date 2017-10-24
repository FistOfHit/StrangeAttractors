% Function to bring together all measure of dimension into one place on one
% 3-D attractor for your veiwing pleasure and convinience.

% Inputs are the parameters of the attractor and no outputs. 

function [] = FractalDimension(a, b, c, d1, d2, irrationalFreq)

    [position, speed, tauVector] = PoincareByParameters(a, b, c, d1, d2, irrationalFreq, 10^6);

    lyapunovDim = LyapunovDimension(a, b, c, d1, d2, irrationalFreq);
    
    correlationDim = CorrelationDimension(a, b, c, d1, d2, irrationalFreq, 3);
    
    sphereCountDim = SphereCountDimension(a, b, c, d1, d2, irrationalFreq, 3);
    
    numPoints = length(position);
    plot3(position(0.00002 * numPoints : numPoints), speed(0.00002 * numPoints : numPoints), ...
        tauVector(0.00002 * numPoints : numPoints), '+', 'MarkerSize', 0.1);
    
    xlabel('x');
    ylabel('dx/dt');
    zlabel('tau');
    
    title(['L - Dim = ' num2str(lyapunovDim) ', C - Dim = '...
        num2str(correlationDim) ', S - Dim = ' num2str(sphereCountDim) '.']);

end