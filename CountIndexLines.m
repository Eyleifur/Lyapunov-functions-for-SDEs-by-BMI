% Fall sem gefur fj�lda samsetninga af n t�lum sem summast upp � 4 �n
% endurtekninga
% Inntak: n - fj�ldi s�ta
function count = CountIndexLines(n)
    if (n < 2)  % Inntak skr�ti� ef �etta gildir
        count = -1;
        return
    end
    count = 0;
    
    % L�nur sem innihalda 4, (3, 1) e�a (2, 2)
    count = count + n;                  % Fj�ldi lei�a til a� setja 4 � n s�ti
    count = count + n * (n - 1);        % Get sett 3 � eitt af n s�tum, s��an 1 � hin n-1 s�tin
    count = count + n * (n - 1) / 2;    % Einn 2 � n s�ti, hinn 2 � n - 1 s�ti. Deili me� tveimur vegna tv�talninga
    
    % L�nur sem innihalda (2, 1, 1)
    if (n >= 3)
        count = count + n * (n - 1) * (n - 2) / 2;  % Set 2 � eitt af n s�tum, annan 1 � hin n-1 og hinn 1 � n-2. Deili vegna tv�talningar
    end
    
    % L�nur sem innihalda (1, 1, 1, 1)
    if (n >= 4)
        count = count + factorial(n) / (factorial(n - 4) * factorial(4));   % n velur 4
    end
end
        