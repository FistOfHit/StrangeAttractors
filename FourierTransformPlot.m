% Function to plot fourier transformed data against the frequencies to see
% peaks, takes inputs required by FourierTransformByParameters

function [] = FourierTransformPlot(a, b, c, d1, d2, irrationalFreq, initConditions)

    % getting the FT'd data
    freqSpectrum = FourierTransformByParameters(a, b, c, d1, d2, irrationalFreq, initConditions);
   
    % Seeting up the X-axis - frequency
    numPoints   = length(freqSpectrum);
    radsPerSec  = 100;
    frequencies = linspace(-numPoints / 2, numPoints / 2, numPoints) * (radsPerSec / numPoints);
    
    plot(frequencies, freqSpectrum);
    xlabel('Frequency(Hz)');
    ylabel('Amplitude');
    title(['plot at d2 = ' num2str(d2) '.']);
    
% can restrict axis to make it clearer if needed
%     axis([1, 10, 0, 0.12]);

end