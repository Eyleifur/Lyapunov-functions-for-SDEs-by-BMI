% Skilar v�si sem v�sar � �ann sta� � M fylkinu �ar sem vi� setjum Q fylki�
% -1
function M = Qstart(n)
    M = 0.5*n*(2*n^2 - 3*n + 3);
    
    if (n >= 4)
        M = M + 2*nchoosek(n,4);
    end
end