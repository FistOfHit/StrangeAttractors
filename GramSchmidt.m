% Modified Gram-Schmidt algorithm to find new orthonormal basis from old
% one (oldBasis -> new Basis). V is the the 3D matrix we use to store
% the vectors being processed. 

function [newBasis] = GramSchmidt (oldBasis)
    
    [numRows, numCols] = size(oldBasis);
    storageMatrix = zeros(numRows, numCols, numCols);
    newBasis = zeros(numRows, numCols);
        
    % Intialising Storage matrix with old basis
    for i = 1 : numCols
        storageMatrix(:, 1, i) = oldBasis(:, i);
    end
    
    % first orthonormal vector
    newBasis(:, 1) = storageMatrix(:, 1, 1) / norm(storageMatrix(:, 1, 1));
    
    % subtracting projection of new basis vectors from the latest addition to new basis
    for i = 2 : numCols
        
        for j = i : numCols
            storageMatrix(:, i, j) = storageMatrix(:, i - 1, j) - ...
                dot(storageMatrix(:, i - 1, j), newBasis(:, i - 1)) * ...
                newBasis(:, i - 1);
        end
        
        newBasis(:, i) = storageMatrix(:, i, i) / (norm(storageMatrix(:, i, i)));
        
    end
    
end