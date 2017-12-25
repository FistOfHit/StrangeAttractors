function[fractal_dim] = Boxcount_dim(a, b, c, force1, force2, irr_freq, num_cycles)
% Find box counting dimension of fractal-like attractor generated with given parameters.
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
% This measure of fractal dimension is calculated by surrounding the fractal object with a grid or 
% mesh of "boxes" of a certain size in 3D and counting the number of boxes that contain a point in 
% the object in them. The change in the number of boxes required for this as the size of the boxes 
% changes is then observed, and has been proven mathemtically to vary as a power, which can be found
% using a log-log plot. The real genius here comes from the Box_count function developed by Dr.
% Ramsden, see its documentation for more information.

    % Generating attractor with parameters given, removing transients
    [pos, speed, tau] = Create_attractor(a, b, c, force1, force2, irr_freq, num_cycles);
    pos = pos(0.2 * num_cycles : num_cycles);
    tau = tau(0.2 * num_cycles : num_cycles);
    speed = speed(0.2 * num_cycles : num_cycles);
    % Storing minimum values in each dimension
    lower_bounds = [min(pos), min(speed), min(tau)];
    
    num_tests = 20;
    box_count = zeros(1, num_tests);
    box_sizes = 10.^linspace(0, -1, num_tests);
    % Counting boxes of certain size with Boxcount function
    for i = 1 : numTests
        [numBoxes] = Box_count(lower_bounds, pos, speed, tau, box_sizes(i));
        box_count(i) = numBoxes;
    end
    
    % Finding gradient of line of best fit on Log-Log data
    box_count = log(box_count);
    inverse_sizes = log(1./box_sizes);
    fractal_dim = polyfit(inverse_sizes, box_count, 1);
    
    plot(inverse_sizes, box_count, '+');
    hold on
    plot(inverse_sizes, (fractal_dim(1).*inverse_sizes) + fractal_dim(2));
    hold off
    
end