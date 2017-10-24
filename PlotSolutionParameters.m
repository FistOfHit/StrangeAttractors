function [] = PlotSolutionParameters (a, b, c, d1, d2, irrationalFreq)

    skipCycles = (10^4) - 10;
    skipIterations = 100 * skipCycles;
    numIterations = (10^4) * 100;
    
    [position1, speed1] = DuffingSolutionByParameters(a, b, c, d1, d2, irrationalFreq, [1, 0]);
%     [position2, speed2] = DuffingSolutionByParameters(a, b, c, d1, d2, irrationalFreq, [1.001, 0]);
%     [position3, speed3] = DuffingSolutionByParameters(a, b, c, d1, d2, irrationalFreq, [1.01, 0]);
    
% -----Plotting time and position (to see periodic behaviour)-------------
    t = linspace(0, 10^4 * 2 * pi, numIterations);
    
    plot(t(skipIterations : numIterations), position1(skipIterations : numIterations));
    xlabel('t');
    ylabel('x');
%     
%     hold on
%     plot(t(skipIterations : numIterations), position2(skipIterations : numIterations));
%     xlabel('t');
%     ylabel('x');
%     
%     plot(t(skipIterations : numIterations), position3(skipIterations : numIterations));
%     xlabel('t');
%     ylabel('x');
%     
%     title(['plot at d2 = ' num2str(d2) '.']);
%     hold off
% ------------------------------------------------------------------------


%-----Plotting position and second derivative (to see attractors)---------
%     plot(position1(skipIterations : numIterations), speed1(skipIterations : numIterations));
%     xlabel('x');
%     ylabel('dx/dt');
%     % end point as red cross to see where it ends up amdist chaos
%     hold on
%     plot(position1(numIterations), speed1(numIterations), 'r+');
%     
%     plot(position2(skipIterations : numIterations), speed2(skipIterations : numIterations));
%     xlabel('x');
%     ylabel('dx/dt');
%     % end point as red cross to see where it ends up amdist chaos
%     plot(position2(numIterations), speed2(numIterations), 'r+');
%     
%     plot(position3(skipIterations : numIterations), speed3(skipIterations : numIterations));
%     xlabel('x');
%     ylabel('dx/dt');
%     % end point as red cross to see where it ends up amdist chaos
%     plot(position3(numIterations), speed3(numIterations), 'r+');
%     hold off
%-------------------------------------------------------------------------
    
end