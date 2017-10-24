% Function to calculate the frequency spectrum of the solution wave, and
% represents it as a double sided spectrum. only performs the transform on
% the last 10*5 values of position to observe frequencies present in long
% term solution only, aswell as a being a convinience for performance and
% calculations

function [freqSpectrum] = FourierTransformByParameters(a, b, c, d1, d2, irrationalFreq, initConditions)

%     test data
%     t = linspace(0, 20000*pi, 10^6);
%     pos = sin(t) + 0.5*sin(4*t) + 0.3*sin(9*t);

    % Initialising values needed
    [position, ~] = DuffingSolutionByParameters(a, b, c, d1, d2, irrationalFreq, initConditions);
    numPoints = length(position);
    
    % calculating the DFT, shifting it to show a positive and negative spectrum
    freqSpectrum = fft(position(numPoints - (10^5) : numPoints));
    freqSpectrum = fftshift(freqSpectrum);
    freqSpectrum = abs((200 * freqSpectrum) / numPoints);

end