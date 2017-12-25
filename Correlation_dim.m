function[fractal_dim] = Correlation_dim(a, b, c, force1, force2, irr_freq, num_cycles)
% Find correlation dimension of fractal-like attractor generated with given parameters
%  
% Inputs:
% a = -1*Damping strength
% b = Linear stiffness
% c = non-linearity in restoring force
% force1 = amplitude of force in rational frequency driving force (2*pi Hz)
% force2 = amplitude of force in irrational frequecy driving force (2*pi*irr_freq) Hz)
% irr_freq = value of irrational driving frequency 
% num_cycles = number of driving cycles used to generate attractor
% 
% Outputs: 
% fractal_dim = fractal dimension calculated via box counting algorithm
% 
% Figures:
% Log-log plot of box count against box size, as well as fitted linear model
%
% This measure of fractal dimension is found by essnetially measuring the probability of points
% being "near" any given point, or in other words, the probability of the attractor having a certain
% density. This is done by looking at each point in the attractor and determining how many other
% points lie within a certain radius of it. We then observe how this quantity changes as the chosen
% radius changes, and mathematically its proven that it varies as a power and so a log-log plot is
% used to extract that as a linear gradient. For more information on the correlation dimension,
% check out the final report in this repository. 
%
% Note: As discussed in the report, the correlation dimension converges VERY slowly, hence any
% accurate result can only be obtained by using 10^4 points in the attractor (num_cycles) or more.
% However, as its an O(N^2) algorithm, good luck with anything more than 10^6. 

    % number of iterations for varying radii
    num_tests = 25;
    
    % Generating data and removing transient effects
    [pos, speed, tau] = Create_attractor(a, b, c, force1, force2, irr_freq, num_cycles);
    
    % Skipping ALL transients in further calculations
    pos = pos((0.2*num_cycles) : num_cycles);
    speed = speed((0.2*num_cycles) : num_cycles);
    tau = tau((0.2*num_cycles) : num_cycles);
    num_points = length(pos);
    
    corr_sums = zeros(1, num_tests);
    radii = 10.^linspace(-1.5, -2.5, num_tests);

    for i = 1 : num_tests

        radius = radii(i);
        num_in_radius = 0;

        for j = 1 : num_points

            % Looking at a certain point in attractor
            p1 = pos(j);
            s1 = speed(j);
            t1 = tau(j);

            % Setting up the rest of the points as a vector...
            p2 = pos(j + 1 : num_points);
            s2 = speed(j + 1 : num_points);
            t2 = tau(j + 1 : num_points);

            % ...So it saves a lot of time by using inbuilt element-wise operators.
            distances = (p1 - p2).^2 + (s1 - s2).^2 + (t1 - t2).^2;

            % Array of bits after a condition to mimic heaviside step function
            is_in_radius  = (radius^2 - distances) > 0;
            num_in_radius = num_in_radius + sum(is_in_radius);

        end 

        % Correlation sum, approximation to correlation intgral.
        corr_sums(i) = (num_in_radius) / (num_points*(num_points - 1));

    end
   
    % setting up recorded values for linear regression
    corr_sums = log(corr_sums);
    radii = log(radii);
    % fractal dimension as gradient via basic linear regression
    fractal_dim = polyfit(radii, corr_sums, 1);
    
    % Plotting attractor and data with linear fit
    figure1 = figure;

    subplot(1, 2, 1);
    plot3(pos, speed, tauVector, '.', 'MarkerSize', 1);
    title('Plot of attractor');
    xlabel('x');
    ylabel('dx/dt');
    zlabel('tau');

    subplot(1, 2, 2);
    plot(radii, corrSums);
    hold on
    plot(radii, fractalDim(1).*radii + fractalDim(2), '+');
    hold off
    title('Number of local points against chekcing distance');
    xlabel('log(radius)');
    ylabel('log(number of points)');
    
    suptitle(['Correlation dimension = ' num2str(fractalDim(1)) '.']);
    
    fractal_dim = fractal_dim(1);
    
end