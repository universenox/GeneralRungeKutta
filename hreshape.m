function[sol] = hreshape(vec, rows, cols)
    % reshapes vector into rows x cols matrix, by row first
    sol = zeros(rows,cols);
    itr = 1;
    
    for i = 1:rows
        for j = 1:cols
            sol(i,j) = vec(itr);
            itr = itr + 1;
        end
    end
end