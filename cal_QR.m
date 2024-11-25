function [Q, R] = cal_QR(eigenS, eigenS_multi)
%%%% Construct Q, E from eigenvalues %%%%

if size(eigenS, 2) ~= size(eigenS_multi, 2)
    error("Eigenvalues and multiplicities should be same length.")
end
n = size(eigenS, 2);
dim_of_blocks = zeros(n);

for i = 1:n
   if imag(eigenS(i)) == 0
       dim_of_blocks(i) = eigenS_multi(i);
   else
       dim_of_blocks(i) = eigenS_multi(i) * 2;
   end
end
v = sum(dim_of_blocks, "all");
S = zeros(v, v);
L = zeros(v, 1);

% Construct Jordan Blocks for S and L
start_idx = 1; 
for i = 1:n
    for j = 1: eigenS_multi(i)
        if j == 1
            L(start_idx) = 1;
        end
        
        if imag(eigenS(i)) == 0
            S(start_idx+j-1, start_idx+j-1) = eigenS(i);
            if j ~= eigenS_multi(i)
                S(start_idx+j-1, start_idx+j) = 1; 
            end
        else
            from_idx = start_idx+ 2*j -2;
            num = imag(eigenS(i));
            S(from_idx: from_idx+1, from_idx: from_idx+1) = [0, num; -num, 0];
            L(from_idx) = 1;
            if j ~= eigenS_multi(i)
                S(from_idx: from_idx+1, from_idx + 2: from_idx+1 + 2) = eye(2);
            end
        end    
    end
    start_idx = start_idx+dim_of_blocks(i);
end

Q = S;
R = L;
end