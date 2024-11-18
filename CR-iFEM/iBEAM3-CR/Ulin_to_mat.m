function matU = Ulin_to_mat(Ulin,width)
    matU = [];
    c = width;
    r = size(Ulin,1)/width;
    for i = 1:r
        for j = 1:c
            k = (i-1)*c+j;
            matU(i,j) = Ulin(k);
        end
    end
end