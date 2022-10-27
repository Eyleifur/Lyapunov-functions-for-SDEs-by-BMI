% Skilar vísi sem vísar á þann stað í M fylkinu þar sem við setjum Q fylkið
% -1
function M = Qstart(n)
    M = 0.5*n*(2*n^2 - 3*n + 3);
    
    if (n >= 4)
        M = M + 2*nchoosek(n,4);
    end
end