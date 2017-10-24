% Function to plot the solutions of the system, taking the inputs required
% for DuffingSolution as well as the number of cycles you want to skip.
% uncomment/comment the code you want to use/not use.

function [] = PlotSolution (stepsPerCycle, skipCycles, numCycles, initConditions)
    
    numIterations  = stepsPerCycle * numCycles;
    skipIterations = stepsPerCycle * skipCycles;
    
    [position, speed] = DuffingSolution(stepsPerCycle, numCycles, initConditions);
    
% -----Plotting time and position (to see periodic behaviour)-------------
    t = linspace(0, numCycles * 2 * pi, numIterations);
    plot(t(skipIterations + 1 : numIterations) / (2*pi), position(skipIterations + 1 : numIterations));
    xlabel('t (scaled by 1/2pi)');
    ylabel('x');
% ------------------------------------------------------------------------


%-----Plotting position and second derivative (to see attractors)---------
%     plot(position(skipIterations + 1 : numIterations), speed(skipIterations + 1 : numIterations));
%     xlabel('x');
%     ylabel('dx/dt');
%     % end point as red cross to see where it ends up amdist chaos
% %     hold on
% %     plot(position(numIterations), speed(numIterations), 'r+');
% %     hold off
%-------------------------------------------------------------------------
    
end