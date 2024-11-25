clear; clc; 
close all;

load Lin_WFModel_20WTG.mat
Lin_InterPoints; % interpolation points (S, L)

A = linsys.A; B = linsys.B; C = linsys.C;
sys = ss(A, B, C, []);

%% ROM
Lin_SignalGenerator;
Lin_Sylvester2Moments;

G = pinv(YPI_syl)*YB_syl;
sysr = ss(S-G*L, G, CPI_syl, 0);

stability = isstable(sysr);

