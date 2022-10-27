% Fall skilar fylki af stærð n með 1 í sætum (i,j) og (j,i) og 0 annars.
% Þar sem E er nxn einingafylki
function Ev = Eij(E,i,j)
    if i == j
        Ev = E(:,i)*E(j,:);
    else
        Ev = E(:,i)*E(j,:) + E(:,j)*E(i,:);
    end
end