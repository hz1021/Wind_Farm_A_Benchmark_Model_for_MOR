%% SignalGenerator 
% Interpolation points
s = linspace(0, 6.1, 10);
upsilon = s + 0.1;
eigenS = 10.^s.*1i;
eigenQ = 10.^upsilon.*1i;

eigenS_multiplicity = ones(1, size(eigenS, 2));
eigenQ_multiplicity = ones(1, size(eigenQ, 2));
