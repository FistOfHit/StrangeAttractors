function[] = DuffingRK4_PlotAnimation1(n, tmax, ic1, ic2)
    
%-----Range of initial position--------------------------------------------
    ic1Range = linspace(ic1(1), ic1(2), 100);
    
%-----Loop to plot animation-----------------------------------------------    
    for i=1:(length(ic1Range) - 1)
        [pos, speed] = DuffingRK4(n, tmax, [ic1Range(i), ic2]);
    
%-----Plot of position and second derivative (changing over time)----------        
        plot(pos, speed);
        xlabel('x');
        ylabel('dx/dt');
        x0 = ic1Range(i);
        
        
        hold on
        plot(pos(n-1), speed(n-1), 'r+');
        title(['x(0) = ' num2str(x0) '.']);
        hold off
        
%-----pause in real time for animation-------------------------------------        
        pause(0.01);
        
    end
    
    text(0,0,'complete');
    
end