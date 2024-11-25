function [vector] = RBFbasisnD(x, x_c, sigma)
%POLYBASIS Summary of this function goes here
%   Detailed explanation goes here
n = size(x, 1);  % [x1; x2; x3]
n_basis = size(x_c, 1); 
vector = zeros(n, n_basis);
% vector = zeros(n, n_basis + 1);
for j  = 1:n
    for i = 1:n_basis
       vector(j, i) = exp(-norm(x(j, :) - x_c(i, :), 2)^2./(2*sigma(i)^2));
    end
end
% vector(:, end) = ones(n, 1);
end
