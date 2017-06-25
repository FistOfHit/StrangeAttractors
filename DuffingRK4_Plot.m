function[] = DuffingRK4_Plot(n, tmax, ic)
    
%-----Initialising arrays--------------------------------------------------
     time = linspace(0, tmax, n);
     [pos, speed] = DuffingRK4(n, tmax, ic);
    
%-----Plotting time and position (see periodic behaviour)------------------
     plot(time, pos);
     xlabel('t');
     ylabel('x');
    
%-----Plotting time and second derivative (dont know what im seeing yet)---
     %plot(time, speed);
     %xlabel('t');
     %ylabel('dx/dt');
    
%-----Plotting position and second derivative (see attractor)--------------
     %plot(pos, speed);
     %xlabel('x');
     %ylabel('dx/dt');
    
    %-----end point as red cross to see where it ends up amdist chaos------
     %hold on
     %plot(pos(n-1), speed(n-1), 'r+');
     %hold off
    
%-----Poincare section plot------------------------------------------------

      
         
         
         
         
     
          
     
     
end
