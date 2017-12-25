function[] = Plot_solutions(a, b, c, force1, force2, irr_freq, init_cond)
% Generate position or phase plots of solutions of duffing oscillator.
%
% Inputs:
% a = -1*Damping strength
% b = Linear stiffness
% c = non-linearity in restoring force
% force1 = amplitude of force in rational frequency driving force (2*pi Hz)
% force2 = amplitude of force in irrational frequecy driving force (2*pi*irr_freq) Hz)
% irr_freq = value of irrational driving frequency 
% init_cond = initial conditions for dynamical system
%
% Outputs:
% None.
%
% Figures:
% Plot of position against time
% Plot of speed against position
%
% The plots created are in the beautiful and standard blue, for the optimum veiwing pleasure of the
% user. Data is generated via the Duffing_solution function, and then all plotted in one figure side
% by side. Due to this being a seperate function to Duffing_solution, some values such as number of
% iterations or time must be calculated here with knowledge of chosen step size etc. from the
% Duffing_solution function. skip_iters is how many iterations of the solution you want to skip from
% the start in order to see the long term chaotic behaviour of the oscillator. Reccomended skip is
% 99.7% or so of the iterations for any reasonable number of cycles (10^4 by default).
    
    % Generating solutions with Duffing_solution function
    [pos, speed] = Duffing_solution(a, b, c, force1, force2, irr_freq, init_cond);
    
    % Number of iterations and number to skip from start in plot (to observe long-term behaviour)
    num_iters = 10^6;
    skip_iters = ceil(9.97*10^5);
    time = linspace(0, 10^4*2*pi, num_iters);
    
    figure1 = figure;
    % Plotting time and position
    subplot(1, 2, 1);
    plot(time(skip_iters + 1 : num_iters), pos(skip_iters + 1 : num_iters));
    title('position against time');
    xlabel('t');
    ylabel('x');

    % Plotting position and speed (Phase plot)
    subplot(1, 2, 2);
    plot(pos(skip_iters + 1 : num_iters), speed(skip_iters + 1 : num_iters));
    title('speed against position (Phase plot)')
    xlabel('x');
    ylabel('dx/dt');
    % end point as red cross to see where it ends up amdist chaos
    hold on
    plot(pos(num_iters), speed(num_iters), 'r+');
    hold off
    
    suptitle('Plots of solutions of QP-driven duffing oscillator');
    
end