function mM = BMISolution(F, G, p0, eps, seed)
    n = length(F);
    % Lesa inn fylki
    numGs = size(G,3);
    if numGs == 1
        load([num2str(n), 'x', num2str(n), 'BMI System.mat'], 'A', 'B', 'C', 'Peqz', 'S');
    else
        load([num2str(numGs), 'G ', num2str(n), 'x', num2str(n), 'BMI System.mat'], 'A', 'B', 'C', 'Peqz', 'S');
    end
    
    syms p epsi
    
    F1 = sym('F%d', [n,n]); %[n,n]
    G1 = sym('G%d', [n,n,numGs]); %[n,n]
    %subs(A(:,:,:,:), {F1, G1, p, epsi}, {F, G, p0, eps})
    A = subs(A, F1, F);                 % Symbolic to numeric
    A = subs(A, G1, G);
    A = subs(A, [p, epsi], [p0, eps]);
    A = double(A);
    B = double(B);
    C = subs(C, epsi, eps);
    C = double(C);
    
    maxIter = 100;                % maximum number of iterations
    eta     = 50;                 % parameter which controls importance of penalty (for a problem of this size, a number between 50-100 works well)
    
    rng(seed);
    Msize = BMISize(n);
    Bsize = size(B,3);
    m = n*(n+1) / 2;
    q0 = rand(Bsize, 1); %rand(Bsize, 1) held ég              % random starting point
    c = rand(Bsize, 1);  %rand(Bsize, 1)             % random linear objective function specified by c
    
    
    for i = 1 : maxIter 

        cvx_begin
        cvx_solver sdpt3
        cvx_quiet true
        cvx_precision best
        variable q(Bsize, 1)               
        variable Q(Bsize, Bsize) symmetric      
        
        M = zeros(Msize);
        for j = 1:m
            for k = j:m
                M = M + Q(j, k)*A(:,:,j,k);
            end
        end
        
        for j = 1:Bsize
            M = M + q(j)*B(:,:,j);
        end
        
        M = M + C;
        %M = Q(1, 1)*A11+Q(1, 2)*A12+Q(1, 3)*A13+Q(2, 2)*A22+Q(2, 3)*A23+Q(3, 3)*A33+q(1)*B1+q(2)*B2+q(3)*B3+q(4)*B4+q(5)*B5+q(6)*B6+C;

        minimize(c'*q + eta*(trace(Q) - 2*q0'*q)) 

        subject to

        M == semidefinite(Msize);
        [Q, q; q',1] == semidefinite( Bsize + 1 );

        cvx_end

        fprintf('Infeasibility: %f \n', trace(Q-q*q'));     
        q0=q;  
    end
    
    format long
    format compact
    
    %mM=q(1)^2*A11+q(1)*q(2)*A12+q(1)*q(3)*A13+q(2)^2*A22+q(2)*q(3)*A23+q(3)^2*A33+q(1)*B1+q(2)*B2+q(3)*B3+q(4)*B4+q(5)*B5+q(6)*B6+0.5*C;
    mM = zeros(Msize);
    for j = 1:m
        for k = j:m
            mM = mM + q(j)*q(k)*A(:,:,j,k);
        end
    end
    for j = 1:Bsize
        mM = mM + q(j)*B(:,:,j);
    end
    mM = mM + C;
    
    Q0 = Qstart(n);
    mQ = mM((Q0+1):(end-1), (Q0+1):(end-1));
    Q = mQ
    eigQ = eig(mQ);
    mP = mM(1:m,1:m);
    eigP = eig(mP)';
    c = mM(Msize, Msize)
    if (sum(~(eigQ >= 0)) + sum(~(eigP >= 0)) == 0)
        fprintf('Solution found\n');
        fprintf('Eigenvalues of Q:\n');
        eigQ
        fprintf('Eigeinvalues of P:\n');
        eigP
    else
        fprintf('No solution found.\n');
        eigQ
        eigP
    end
end