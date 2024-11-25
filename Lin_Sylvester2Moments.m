%   Input: S, A, M, C
%   Return: Upsilon * B

% Compute PI
PI_syl = sylvester(A, -S, -B*L);
CPI_syl = C * PI_syl;

% Compute Y
Y_syl = sylvester(Q, -A, R*C);
YB_syl = Y_syl * B;

% Compute YPI
YPI_syl = Y_syl * PI_syl;
