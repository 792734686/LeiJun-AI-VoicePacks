function linspaceU = Umat_to_lin(Umat)
linspaceU = [];
    [r,c] = size(Umat);
    for i = 1:r
        for j = 1:c
            k = c*(i-1)+j;
            linspaceU(k,1) = Umat(i,j);
        end
    end
end