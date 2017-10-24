function [fractalDim] = BoxCountDimension(a, b, c, d1, d2, irrationalFreq, numCycles)

    [position, speed, tauVector] = PoincareByParameters(a, b, c, d1, d2, irrationalFreq, numCycles);
    position  = position(0.2 * numCycles : numCycles);
    speed     = speed(0.2 * numCycles : numCycles);
    tauVector = tauVector(0.2 * numCycles : numCycles);
    
    lowerBounds = [min(position), min(speed), min(tauVector)];

    numTests = 20;
    
    boxCount = zeros(1, numTests);
    boxSizes = 10 .^ linspace(0, -1, numTests);
    
    for i = 1 : numTests
        [numBoxes] = ProcessData(lowerBounds, position, speed, tauVector, boxSizes(i));
        boxCount(i) = numBoxes;
    end

    boxCount = log(boxCount);
    invBoxSizes = log(1 ./ boxSizes);
    
    fractalDim = polyfit(invBoxSizes, boxCount, 1);
%     fractalDim = fractalDim(1);
    
    plot(invBoxSizes, boxCount, '+');
    hold on
    plot(invBoxSizes, (fractalDim(1) .* invBoxSizes) + fractalDim(2));
    hold off
    
end