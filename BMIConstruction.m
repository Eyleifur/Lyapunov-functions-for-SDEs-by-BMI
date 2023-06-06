function [A, B, C, S, Peqz] = BMIConstruction(n, NumGs)
    % Skilgreinum breytur

    syms p c epsi;
    Q = sym('Q', [n,n]);
    Q = triu(Q, 0) + triu(Q, 1).'; % Samhverft fylki

    F = sym('F%d', [n,n]);
    
    % Ef ekkert er sett inn fyrir NumGs þá gerum við ráð fyrir einu G-fylki
    if (nargin == 1)
        G = sym('G%d', [n, n, 1]);
    else
        G = sym('G%d', [n, n, NumGs]);
    end
    
    xv = sym('x%d', [n, 1]);

    Zf = xv * xv.';

    Z = tril(Zf, 0).';

    m = n*(n+1) / 2;
    Z = Z(find(Z));
    
    GQG = 0;
    xQG = 0;
    for i = 1:NumGs
        GQG = GQG + G(:,:,i).'*Q*G(:,:,i);
        xQG = xQG + (xv.'*(Q*G(:,:,i) + G(:,:,i).'*Q)*xv)^2;
    end

    Pc = -xv.'*(F.'*Q + Q*F + GQG)*(xv*xv.')* Q*xv + (2 - p)/4 * xQG - c*(xv.'*xv)^2;

    P = sym('P', [m,m]);
    P = triu(P, 0) + triu(P, 1).';
    PZ = Z.' * P * Z;


    nn = CountIndexLines(n);
    I = zeros(nn, n); % Index fylki sem geymir hversu oft diffrað er m.t.t hvers staks. dx(1)^I(1)...dx(n)^I(n) þar sem sum(I) = 4

    count = 1; % Vísir á röð í fylki I

    % Bý til allar samsetningar af tölum (i_1, ... i_n) þannig að
    % sum(i_1,...,i_n) = 4. Þetta verða vísarnir í diffruninni.
    % Athugum einnig hvenær tvær eða fleiri eins tölur koma fyrir í sömu
    % línunni og geymum í fylkinu dep. Þurfum það fylki seinna
    for i=1:n
        for j=i:n
            for k=j:n
                for l=k:n
                    I(count, i) = I(count, i) + 1;
                    I(count, j) = I(count, j) + 1;
                    I(count, k) = I(count, k) + 1;
                    I(count, l) = I(count, l) + 1;
                    count = count + 1;
                end
            end
        end
    end
    I;

    % Finnum þær línur sem innihalda háðar P breytur
    dep = zeros([nn 1]);
    numDep = 0;
    for i = 1:nn
        b = histc(I(i,:), [0 1 2 3 4]);
        if (b(2) > 1) ||(b(3) > 1) || (b(4) > 1) || (b(5) > 1)
            dep(i) = 1;
            numDep = numDep + 1;
        end
    end

    Peqz = sym('a', [1 nn]);    % Vigur sem inniheldur vinstri hlið sambandanna, þ.e. stuðlar PZ
    f = sym('a', [1 nn]);       % Vigur sem inniheldur hægri hlið sambandanna, þ.e. stuðlar Pc

    fact = [1, 1, 2, 6, 24]; % Forreiknum i! fyrir i = 0,...,4

    % Diffrum PZ fjórum sinnum m.t.t. til allra samsetninga af breytum í vigrinum xv.
    % Ytri for-lykkjan fer í gegnum allar línur í index fylkinu I.
    % Innri for-lykkjan fer í gegnum breyturnar í xv.
    % Stakið (i,j) í I segir að í i-ta sambandi þá er PZ diffrað I(i,j) sinnum
    % með tilliti til breytunnar xv(j).

    % Deilum í hverri ítrun með fact(I(i,j) + 1) til að losna við fastana sem fylgir
    % diffruninni.
    % Ef í i-tu línu í I er að finna ás þá
    % þarf að deila aukalega með tveimur

    % Diffrum Pc í sömu for-lykkju og geymum í vigrinum f
    for i = 1:nn
        g = PZ; % Dummy breyta svo PZ breytist ekki
        h = Pc;
        ace = 0; %  Athugum hvort ás sé að finna í i-tu línu
        for j = 1:n
            g = diff(g, xv(j), I(i,j)) / fact(I(i,j) + 1);
            h = diff(h, xv(j), I(i,j)) / fact(I(i,j) + 1);
            if (I(i, j) == 1)
                ace = 1;
            end
        end
        Peqz(i) = g / (1 + ace); % Geymum niðurstöður. Ef ás er í línu þá er ace = 1 og við deilum með tveimur. Annars er ace = 0 og við deilum með 1
        f(i) = h / (1 + ace);
    end
    Peqz.';
    f.';

    % Diffrum nú hvert stak í f til að fá stuðlana við breyturnar Q_iQ_j, i,j =
    % 1,...,n.
    % Fyrst setjum við stökin úr efra hornalínufylki Q í vigurinn q til að
    % einfalda vísanir

    q = Q(triu(true(size(Q)), 0));

    m = n * (n + 1) / 2;    % stærð vigurs q

    % Í Pc þá koma breyturnar Q_i bara fyrir sem annars stigs liðir svo
    % nægjanlegt er að diffra tvisvar hvert stak í f.
    % Geymum niðurstöður í fylkinu S.
    % Athugum að c er einnig breyta svo við þurfum einnig að diffra m.t.t. c.
    % Bætum því við einum dálk í S
    mm = m * (m + 1) / 2;
    S = sym('S', [nn (mm + 1)]);
    m;
    for i = 1:nn
        g = f(i);
        count = 1;
        for j = 1:m
            for k = j:m
                S(i, count) = diff(diff(g, q(j)), q(k)) / (1 + (j == k));   % Ef Q_j^2 er diffrað tvisvar kemur fasti sem þarf að deila í burt
                count = count + 1;
            end
        end
        S(i, count) = diff(f(i), c, 1);
    end
    S;

    % Vídd fylkisins í BMI ójöfnunni
    Msize = BMISize(n);

    % Búum til E fylkin sem við notum til að framleiða fylkin A_ij, B_i og C í
    % fylkjaójöfnunni.

    % Gerum þau út frá fylkinu
    E = eye(Msize);

    % Setjum nú upp fylkjaójöfnuna

    % Byrjum á A-fylkjunum
    % Hvert A-fylki er Msize x Msize af stærð
    % Þurfum jafnmörg A-fylki og fjöldi Q_iQ_j breyta
    % Geymum fylkin í fylki af fylkjum

    A = sym('A', [Msize, Msize, m, m]);
    A(:,:,:,:) = 0;                     % Setjum allar breyturnar í fylkjunum jafnt 0

    for i = 1:m
        for j = i:m
            p = QQVariableIndex(m,i,j);
            count = 0;                          % Telur hve margar línur með háðum breytum hafa komið
            for k = 1:nn
                if (dep(k) == 0)                % Ef lína með óháðum breytum þá stuðull settur í staðsetningu P_ij
                    [row, col] = getPIndex(P, Peqz, k);
                    if row == col
                        A(:,:,i,j) = A(:,:,i,j) + (E(:,row)*E(row,:)) * S(k,p);
                    else
                        A(:,:,i,j) = A(:,:,i,j) + (E(:,row)*E(col,:) + E(:,col)*E(row,:)) * S(k,p);
                    end
                else                            % Annars stuðull settur á hornalínu M. Bæði + og -
                    v = m + 2*count + 1;          % Vísir á það hvert stuðlar háðra breyta eru settir í fylkið. v = size(P) + 2*(fjöldi háðra lína sem hafa komið) + 1
                    A(:,:,i,j) = A(:,:,i,j) + (E(:,v)*E(v,:)) * S(k,p);
                    A(:,:,i,j) = A(:,:,i,j) - (E(:,v+1)*E(v+1,:)) * S(k,p);
                    count = count + 1;
                end
            end
        end
    end

    % Finnum fjölda háðra P-breyta í hverri línu
    depInLine = zeros([nn 1]);
    for i = 1:nn
        reg = regexp(char(Peqz(i)), 'P', 'Match');
        depInLine(i) = length(reg);
    end

    % Byggjum upp B-fylkin
    % Fjöldi B-fylkja = Fjöldi ólíkra staka í Q + Fjöldi óháðra breyta + c
    % = m + 2*(Fjöldi lína af gerðinni (2,2)) + 2*(Fjöldi lína af gerð (2,1,1))
    % + 3*(Fjöldi lína af gerð (1,1,1,1)) + 1

    Bsize = NumBMatrices(n);

    B = sym('B', [Msize Msize Bsize]);
    B(:,:,:) = 0;

    Q0 = Qstart(n);
    numBm = 1; % Vísir á B-fylki til þæginda
    for i = 1:n
        for j = i:n
            p = QQVariableIndex(n,i,j);
            B(:,:,p) = Eij(E, Q0+i, Q0+j);
            numBm = numBm + 1;
        end
    end

    % Þurfum að fiska út indexa háðu breytanna í Peqz
    % Notum regex
    % Gildin í Peqz eru á forminu
    % ()*P()_() + P()_()
    % P()_() + P()_()
    % ()*P()_() + ()*P()_()
    % P()_() + P()_() + P()_()
    % Þar sem () táknar tölu.

    count = 0;
    for i = 1:nn
        if (dep(i) == 1)
            v = m + 2*count + 1;            % Vísir á það hvert stuðlar háðra breyta eru settir í fylkið
            str = char(Peqz(i));
            L = regexp(str, '(?<constant>\d*)\**P(?<first>\d*)_(?<second>\d*)', 'names');
            Num = zeros(3,1);
            for jj = 1:length(L)
                if isempty(L(jj).constant)
                    Num(1) = 1;
                else
                    Num(1) = str2double(L(jj).constant);
                end
                Num(2) = str2double(L(jj).first);
                Num(3) = str2double(L(jj).second);
                
                B(:,:,numBm) = Eij(E,Num(2),Num(3)) + Num(1)*(Eij(E, v+1, v+1) - Eij(E, v, v));
                numBm = numBm + 1;
            end

            count = count + 1;
        end
    end

    % Síðasta B-fylkið. Þetta svarar til breytunnar c
    count = 0;
    for i = 1:nn
        if (dep(i) == 0)
            [row, col] = getPIndex(P, Peqz, i);
            B(:,:,numBm) = B(:,:,numBm) + Eij(E, row, col) * S(i, mm + 1);  % mm + 1-sta stakið er það síðasta. Þar er c
        else
            v = m + 2*count + 1;            % Vísir á það hvert stuðlar háðra breyta eru settir í fylkið
            B(:,:,numBm) = B(:,:,numBm) + Eij(E, v, v) * S(i, mm + 1);
            B(:,:,numBm) = B(:,:,numBm) - Eij(E, v+1, v+1) * S(i, mm + 1);
            count = count + 1;
        end
    end

    % Munum að við geymum breytuna sem svarar til c í síðasta hornalínustaki
    % fylkisins M.
    B(:,:,numBm) = B(:,:,numBm) + Eij(E, Msize, Msize);

    % Loks er það fylki C
    % Hér setjum við inn skilyrðin: c - eps >= 0 og Q - eps*I >= 0.

    C = sym('C', [Msize Msize]);
    C(:,:) = 0;

    for i = (Q0 + 1):Msize
        C(:,:) = C(:,:) + Eij(E, i, i);
    end
    C(:,:) = -epsi*C(:,:);
end