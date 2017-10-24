idfRange = linspace(0.27, 0.3, 20);

a = -0.2; b = 1; c = -1.8; d1 = 0.3;

numCycles = 10^4;

for i = 1 : 1000
    
    [pos1, spd1, tau1] = PoincareByParameters(a, b, c, d1, idfRange(i), ((sqrt(5) - 1) / 2), numCycles);
    exponents = LyapunovExponent(a, b, c, d1, idfRange(i), ((sqrt(5) - 1) / 2));
 
     plot3(tau1((0.2 * numCycles) : numCycles), pos1((0.2 * numCycles) : numCycles),...
        spd1((0.2 * numCycles) : numCycles), '.', 'MarkerSize', 1);
    
    title(['Lexp1 = ' num2str(exponents(1)) ', Lexp2 = ' ...
        num2str(exponents(2)) ', d2 = ' num2str(idfRange(i)) '.']);
    
    pause(0.000001);
    
end