%% Signal Generator
% Compute S, L, Q, and R
[S, L] = cal_QR(eigenS, eigenS_multiplicity);
[Q, R] = cal_QR(eigenQ, eigenQ_multiplicity);
L = L';
nu = size(S, 1);

% Check S and A does not have sharing eigenvalues
[~,D] = eig(A); 
eigenA = diag(D);

for i = 1: size(eigenS, 2)
    for j = 1: size(eigenA, 1)
        if eigenA(j) == eigenS(i) || eigenA(j) == -eigenS(i) || eigenA(j) == eigenQ(i) || eigenA(j) == -eigenQ(i)
            error("S, A should not have sharing eigens")
        end
    end
end
clear i j V D