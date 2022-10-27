% Skilar fjölda B-fylkja í fylkjaójöfnunni
function M = NumBMatrices(n)
    M = 0.5*(2*n^3 - 3*n^2 + 3*n + 2);
    
    if (n >= 4)
        M = M + 3*nchoosek(n, 4);
    end
end