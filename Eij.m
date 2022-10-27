% Fall skilar fylki af st�r� n me� 1 � s�tum (i,j) og (j,i) og 0 annars.
% �ar sem E er nxn einingafylki
function Ev = Eij(E,i,j)
    if i == j
        Ev = E(:,i)*E(j,:);
    else
        Ev = E(:,i)*E(j,:) + E(:,j)*E(i,:);
    end
end