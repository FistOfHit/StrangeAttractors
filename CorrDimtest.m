
numPoints = linspace(3.1, 5, 15);
corrDims = zeros(1, 15);

for i = 1 : length(numPoints)
    corrDims(i) = CorrelationDimension(-0.2, 1, -1.8, 0.3, 0.05, ((sqrt(5) - 1) / 2), floor(10^(numPoints(i))));
end

plot(numPoints, corrDims);
