function [numBoxes] = ProcessData(lowerBounds, position, speed, tauVector, boxSize)

    posMin = lowerBounds(1);
    spdMin = lowerBounds(2);
    tauMin = lowerBounds(3);

    posProcessed = floor((position  - posMin) / boxSize);
    spdProcessed = floor((speed     - spdMin) / boxSize);
    tauProcessed = floor((tauVector - tauMin) / boxSize);

    dataProcessed = transpose([posProcessed; spdProcessed; tauProcessed]);
    uniquePoints = unique(dataProcessed, 'rows');
    numBoxes = length(uniquePoints);

end