clear; clc;
close all

%% Load data
Non_ParameterInitialize;
load("Wt.mat");
load("Yt.mat");

%% Wind speed
ps = 7.5;

%% Approximate using time-domain samples
disp("==> Simulating the data ...");
K = 1; % number of sampling vallues of parameter
N = 40; % number of basis functions
h = 40 * nu; % width of time window 180
sampling = 21000; % interval between two sampling time instant 25
offset = 0 * sampling; 

O = zeros(h*K, 1);
W = zeros(h*K, nu+1); 

p = ps;
wt = Wt;
yt = Yt;
W(1:h, :) = [wt(end-h*sampling+1-offset:sampling:end-offset, :), ones(h, 1) * p]; 
O(1:h, :) = yt(end-h*sampling+1-offset:sampling:end-offset, :); 

disp("done");

W_max = max(W); 
W_min = min(W);
W_between = W_max - W_min; 
% Parameters of basis functions
centers = zeros(N, nu + 1); 
sigmas = ones(N, 1) * 1;
for i_b = 1:N
% centers = linspace(LOW, HIGH, N)';
    for j_b = 1:nu+1
        % centers(i_b, j_b) = rand() * W_between(j_b) + W_min(j_b);
        centers(i_b, j_b) = rand() * 3 * W_between(j_b) + W_min(j_b) - W_between(j_b);
    end 
end
% load("centers_80b0.6912s.mat");
% [W, m, s] = normalizedata(W, 0, 0); 
% W(:, 1) = 0;
U = RBFbasisnD(W, centers, sigmas);
disp("==> Training the approximator ...");
coeffs = U \ O;
% CPI_p_predict_normalized = RBFbasisD(ps_', centers, sigmas) * coeffs;
disp("done");