% Function to calculate Lyapunov Dimension of attractor. Inputs are the
% parameters of the attractor, and output is the fractal dimension, or
% kaplan-yorke dimension. 

% Function works by applying the Kaplan-Yorke conjecture. Lyapunov
% exponents are calculated from the paramters, ordered in terms of size
% (not absolute size) as L1, L2, ..., Ln and them summed, until the sum
% drops below 0. we record how many exponents have been summed as j, and
% the conjecture states that the dimension is: (j - 1) + (sum / |L(j+1)|).

% Quite simple to calculate so more precision can be put into the
% calcuation of the exponents themselves. However, as it is dependent on
% the exponents, we cannot seperate the dimension as we did for
% correlation dimension or sphere counting dimension. Instead, if d2 = 0,
% indicating tau has no effect on the attractor, the dimension is
% automatically dropped by 1, simulating the non-existance of the extra 0
% exponent from tau.

function [fractalDim] = LyapunovDimension(a, b, c, d1, d2, irrationalFreq)

    [position, speed, tauVector] = PoincareByParameters(a, b, c, d1, d2, irrationalFreq, 10^5);
    
    exponents = LyapunovExponent(a, b, c, d1, d2, irrationalFreq);
    orderedExponents = sort(exponents(1 : 3), 'descend');
    disp(orderedExponents)
    
    % finding the sum and the index when the sum becomes negative
    sum = 0;
    for i = 1 : 3
        
        sum = sum + orderedExponents(i);
        
        if sum < 0
            sum = sum - orderedExponents(i);
            break;
        end
        
    end
    
    % by Kaplan - Yorke conjecture
    fractalDim = (i - 1) + (sum / abs(orderedExponents(i)));
    
    numPoints = length(position);
    
    % if d2 = 0, we pretend tau never existed, as it only makes a
    % difference of +- exactly 1 to out result.
    if d2 == 0
        
        fractalDim = fractalDim - 1;
        
        plot(position(0.02 * numPoints : numPoints), speed(0.02 * numPoints : numPoints),...
            '+', 'MarkerSize', 0.1);
        
        xlabel('x');
        ylabel('dx/dt');
        
    else 

        plot3(position(0.02 * numPoints : numPoints), speed(0.02 * numPoints : numPoints),...
            tauVector(0.02 * numPoints : numPoints), '+', 'MarkerSize', 0.1);

        xlabel('x');
        ylabel('dx/dt');
        zlabel('tau');
    
    end
    
    title(['Lyapunov dimension = ' num2str(fractalDim) '.']);
    
end