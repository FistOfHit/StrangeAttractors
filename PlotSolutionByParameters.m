function [] = PlotSolutionByParameters (aRange, bRange, cRange, d1Range, d2Range, IFRange)

    n = 2.5 * (10 ^ 6);
    t = linspace(19980*pi, 20000*pi, 500);
    
    for i = 0:10
        
        % calculating the change in parameters
        x  = i / 10;
        
        a  = (aRange(1)*(1 - x)) + (aRange(2)*x);
        b  = (bRange(1)*(1 - x)) + (bRange(2)*x);
        c  = (cRange(1)*(1 - x)) + (cRange(2)*x);
        d1 = (d1Range(1)*(1 - x)) + (d1Range(2)*x);
        d2 = (d2Range(1)*(1 - x)) + (d2Range(2)*x);
        irrationalFreq = (IFRange(1)*(1 - x)) + (IFRange(2)*x);
        
        [pos, speed] = RK4VectorOutputParameters(a, b, c, d1, d2, irrationalFreq);

        % Plot positioned in center of screen roughly for visibility
        fig1 = figure;
        set(fig1, 'Position', [200, 20, 800, 800])
        
        
        % position againt time
        subplot(2, 2, 1)
        
        plot(t, pos((n - 500) : n));
        xlabel('Time (t)');
        ylabel('Position (x)');
        
        
        % Phase plot
        subplot(2, 2, 2)
        
        plot(pos((n - 500) : n), speed((n - 500) : n));
        
        hold on
        plot(pos(n), speed(n), 'r+');
        hold off
        
        xlabel('Position (x)');
        ylabel('Speed (dx/dt)');
        
            
        % Poincare section
        subplot(2, 2, 3)
        
        [poincareX, poincareY] = RK4PsectionParameters(a, b, c, d, d2);
        
        plot(poincareX, poincareY, 'b.', 'MarkerSize', 3);
        
        xlabel('x');
        ylabel('dx/dt');
        
        
        % Fourier transform to see the underlying dominant frequencies
        subplot(2, 2, 4)
        
        FourierTransformByParameters(a, b, c, d1, d2, 5);
        
        % Title for entire figure
        suptitle(['e = ' num2str(x) ' a = ' num2str(a) ' b = ' num2str(b)...
            ' c = ' num2str(c) ' d = ' num2str(d) ' w = ' num2str(d2)]);
        
        saveas(fig1, ['e = ' num2str(x) '.png']);
        
    end

end