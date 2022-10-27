function [row, col] = getPIndex(P, Peqz, i)
    lin = find(P == Peqz(i),1);        % find finnur alla línulega vísa sem vísa í stakið Peqz(i) í P. Við tökum þann fyrsta (
    n = size(P,1);                     % gefur vídd fylkisins P.
    [col, row] = ind2sub([n n], lin);         % P er samhverft svo við skiptum á col/row þannig að row <= col
end