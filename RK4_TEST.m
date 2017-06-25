function[x] = RK4_TEST(n, tmax, ic1, ic2)

%-----initialising arrays--------------------------------------------------
    t = linspace(0, tmax, n);
    x = zeros(1, length(t));
    y = zeros(1, length(t));
    
%-----Initial conditons and timestep---------------------------------------    
    x(1) = ic1;
    y(1) = ic2;
    h = (tmax) / length(t);
    
%-----dx/dt = F(t,x,y), dy/dt = G(t,x,y)-----------------------------------
    F = @(t,x,y) y;                                 
    G = @(t,x,y) -y;
    
%-----Main calculation loop------------------------------------------------
    for i = 1:(length(t) - 1)
        
        K1 = F(t(i), x(i), y(i));
        L1 = G(t(i), x(i), y(i));
        
        K2 = F(t(i) + (0.5*h), x(i) + (0.5*h*K1), y(i) + (0.5*h*L1));
        L2 = G(t(i) + (0.5*h), x(i) + (0.5*h*K1), y(i) + (0.5*h*L1));
        
        K3 = F(t(i) + (0.5*h), x(i) + (0.5*h*K2), y(i) + (0.5*h*L2));
        L3 = G(t(i) + (0.5*h), x(i) + (0.5*h*K2), y(i) + (0.5*h*L2));
        
        K4 = F(t(i) + h, x(i) + h*K3, y(i) + h*L3);       
        L4 = G(t(i) + h, x(i) + h*K3, y(i) + h*L3);
       
        x(i+1) = x(i) + (1/6)*(K1 + (2*K2) + (2*K3) +K4)*h;
        y(i+1) = y(i) + (1/6)*(L1 + (2*L2) + (2*L3) +L4)*h;
        
    end 
    
end
