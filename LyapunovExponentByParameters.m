% Function to plot the lyapunov exponent as it cahnges with the value of a
% certain parameter or mulitple parameters. can only select one range so
% far but can easily be edited for multiple ranges. Takes 5000 equally
% spaced samples in the range and also plots a line to show where zero is
% and where exactly a system becomes chaotic/non - chaotic.

function [] = LyapunovExponentByParameters(parameterRange)

    precision = 100;
    
    exponentValues = zeros(1, precision);
    parameterValues = linspace(parameterRange(1), parameterRange(2), precision);
    zeroLine = zeros(1, precision);
    
    for i = 1 : precision
        
       lExponents = LyapunovExponent(-0.2, 1, -1.65, 0.3, parameterValues(i), ((sqrt(5) - 1) / 2));
       exponentValues(i) = lExponents(1);
       
    end
    
    % plot the values along with zero-line
    plot(parameterValues, exponentValues, parameterValues, zeroLine);
    
    xlabel('Value of parameter');
    ylabel('Lyapunov exponent of system');
    title('Value of lyapunov exponent as paramter values change');
    
end