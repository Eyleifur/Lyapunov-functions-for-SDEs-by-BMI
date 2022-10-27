function M = BMISize(n)
    M = 0.5*(2*n^3 - 3*n^2 + 5*n + 2);
    
    if (n >= 4)
        M = M + 2*nchoosek(n, 4);
    end
end