function [row, col] = getPIndex(P, Peqz, i)
    lin = find(P == Peqz(i),1);        % find finnur alla l�nulega v�sa sem v�sa � staki� Peqz(i) � P. Vi� t�kum �ann fyrsta (
    n = size(P,1);                     % gefur v�dd fylkisins P.
    [col, row] = ind2sub([n n], lin);         % P er samhverft svo vi� skiptum � col/row �annig a� row <= col
end