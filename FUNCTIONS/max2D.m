function [m,indY,indX] = max2D(A)
    if length(size(A))>2||length(size(A))<2
        error('size of the matrix must be equal to 2!!')
    end
    [mtmp,indY] = max(A);
    [m,indX] = max(mtmp);
    indY = indY(indX);
end