function [chaoticness] = ZeroOneTest(a, b, c, d1, d2, irrationalFreq)
    
    [position, ~, ~] = PoincareByParameters(a, b, c, d1, d2, irrationalFreq, (10^5));
    position = position((9.9 * 10^4) + 1 : 10^5);
%    position = sin(linspace(0, 2000*pi, 1000));

    numPoints = length(position);
    
	numTests = 100;
    cValues = linspace((pi / 5), (pi * (4 / 5)), 100);
    
    cosVar = zeros(1, numPoints);
    sinVar = zeros(1, numPoints);
    
    nCut = numPoints / 10;
    
    meanSquareDisp = zeros(1, nCut);
    
    chaoticness = zeros(1, numTests);
     
    for i = 1 : numTests
        
        c = cValues(i);

        cosVar(1) = position(1) * cos(c);
        sinVar(1) = position(1) * sin(c);

        for j = 2 : numPoints

            cosVar(j) = cosVar(j - 1) + (position(j) * cos(j * c)); 

            sinVar(j) = sinVar(j - 1) + (position(j) * sin(j * c)); 

        end

        for j = 1 : nCut 
			
			cosVarJtoN = cosVar(j : numPoints);
			sinVarJtoN = sinVar(j : numPoints);
            
            cosVar1toJ = cosVar(1 : numPoints - j + 1);
            sinVar1toJ = sinVar(1 : numPoints - j + 1);
			
			displacement = ((cosVarJtoN - cosVar1toJ) .^ 2) + ((sinVarJtoN - sinVar1toJ) .^ 2);
    
            meanSquareDisp(j) = sum(displacement);
            
	    end

        meanSquareDisp = meanSquareDisp / numPoints;
        
        trigInputs = cos(linspace(1, nCut, nCut) * c);
        timeMean   =  sum(position) / numPoints;
        
        oscillatoryError = (1 - trigInputs) * ((timeMean ^ 2) / (1 - cos(c)));
        correctedMSD     = meanSquareDisp - oscillatoryError;

        cutSizes = linspace(1, nCut, nCut);

        correlation = corrcoef(cutSizes, correctedMSD);
        chaoticness(i) = correlation(1, 2);
		
    end
    
    chaoticness = median(chaoticness);
    
end