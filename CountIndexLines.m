% Fall sem gefur fjölda samsetninga af n tölum sem summast upp í 4 án
% endurtekninga
% Inntak: n - fjöldi sæta
function count = CountIndexLines(n)
    if (n < 2)  % Inntak skrítið ef þetta gildir
        count = -1;
        return
    end
    count = 0;
    
    % Línur sem innihalda 4, (3, 1) eða (2, 2)
    count = count + n;                  % Fjöldi leiða til að setja 4 í n sæti
    count = count + n * (n - 1);        % Get sett 3 í eitt af n sætum, síðan 1 í hin n-1 sætin
    count = count + n * (n - 1) / 2;    % Einn 2 í n sæti, hinn 2 í n - 1 sæti. Deili með tveimur vegna tvítalninga
    
    % Línur sem innihalda (2, 1, 1)
    if (n >= 3)
        count = count + n * (n - 1) * (n - 2) / 2;  % Set 2 í eitt af n sætum, annan 1 í hin n-1 og hinn 1 í n-2. Deili vegna tvítalningar
    end
    
    % Línur sem innihalda (1, 1, 1, 1)
    if (n >= 4)
        count = count + factorial(n) / (factorial(n - 4) * factorial(4));   % n velur 4
    end
end
        