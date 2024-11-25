clear; clc;
close all

%% Load data
Non_ParameterInitialize;
load("wt_20WT_78910WS.mat");
load("yt_20WT_78910WS.mat");

%% Hyper-Parameters
LOW = 6;
HIGH = 10;

%% Sample p values
% method 1 - grid sampling
n_samples = 17; 
ps = linspace(LOW, HIGH, n_samples);
% method 2 - Bayesian Regression

%% Approximate using time-domain samples
disp("==> Simulating the data ...");
K = n_samples; % number of sampling vallues of para
N = 40; % number of basis functions
h = 40 * nu; % width of time window 180
sampling = 21000; % interval between two sampling time instant 25
offset = 0 * sampling; 

O = zeros(h*K, 1);
W = zeros(h*K, nu+1); 

for i_p = 1:n_samples
    p = ps(i_p);
    % wt = wt_20WT_400s{i_p};
    % yt = yt_20WT_400s{i_p};
    wt = wt_20WT_78910WS{i_p};
    yt = yt_20WT_78910WS{i_p};
    W(h*(i_p-1)+1:h*i_p, :) = [wt(end-h*sampling+1-offset:sampling:end-offset, :), ones(h, 1) * p]; 
    O(h*(i_p-1)+1:h*i_p, :) = yt(end-h*sampling+1-offset:sampling:end-offset, :); 
end
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

%% Figures 1
% Set wind speed
i_p = 7;
p = ps(i_p);

% Set sampling settings
offset2 = 0 * sampling; 
sampling2 = 50;

% wt = wt_20WT_400s{i_p};
% yt = yt_20WT_400s{i_p};
wt = wt_20WT_78910WS{i_p};
yt = yt_20WT_78910WS{i_p};
wt = wt(1:sampling2:end-offset2, :);
yt = yt(1:sampling2:end-offset2);

% %%%%%%%%
% % Set wind speed
% i_p = 2;
% p = 7.75;
% % Set sampling settings
% offset2 = 0 * sampling; 
% sampling2 = 50;
% 
% wt = wt_FOM_windSpeed{i_p};
% yt = yt_FOM_windSpeed{i_p};
% wt = wt(1:sampling2:end-offset2, :);
% yt = yt(1:sampling2:end-offset2);
% %%%%%%%%

n_t = size(wt, 1);
wp = [wt ones(n_t, 1)*p];
t = linspace(0, 200, n_t);
% [wp, m, s] = normalizedata(wp, m, s);
% wp(:, 1) = 0;
%%%%%%%%%%%%%%%%%%%% Not Interpolation Points %%%%%%%%%%%%%%%%%%%
%%% p = 9.3
% p = 9.3;
% wt = wt_20WT_9_3WS(1:sampling2:end-offset2, :);
% n_t = size(wt, 1);
% wp = [wt ones(n_t, 1)*p];
% yt = yt_20WT_9_3WS(1:sampling2:end-offset2);
% t = t_20WT_9_3WS(1:sampling2:end-offset2);

%%% p = 6.4
% p = 6.4;
% wt = wt_20WT_6_4WS(1:sampling2:end-offset2, :);
% n_t = size(wt, 1);
% wp = [wt ones(n_t, 1)*p];
% yt = yt_20WT_6_4WS(1:sampling2:end-offset2);
% t = t_20WT_6_4WS(1:sampling2:end-offset2);
%%%%%%%%%%%%%%%%%%%%          End             %%%%%%%%%%%%%%%%%%%

kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;

figure;
plot(t, yt);
hold on
plot(t, kappat);
hold off
grid
legend("FOM", "ROM");


% %% Figures 2
% % Set wind speed
% p = 6.4;
% % Set sampling settings
% offset2 = 0 * sampling; 
% tspan = 200; % end time of ROM
% sampling2 = 5;
% 
% out = sim("wind1000_BUS_FOM_ROM");
% t = out.tout;
% t = t(1:sampling2:end-offset2);
% ksit = squeeze(out.Ksait.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p];
% % [wp, m, s] = normalizedata(wp, m, s);
% % wp(:, 1) = 0;
% kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;
% 
% figure;
% % plot(Yt_9_3WS.Time(1:50:end-offset2), Yt_9_3WS.Data(1:50:end-offset2));
% plot(Yt_6_4WS.Time(1:50:end-offset2), Yt_6_4WS.Data(1:50:end-offset2));
% hold on
% plot(t, kappat);
% xlim([80 200]);
% grid on
% legend("FOM", "ROM");
% 
% 
% % %%
% % tspan_set = [40 40 40]; 
% % sampling2 = 50;
% % p_set = [7.5 9.25 6.25]; 
% % 
% % ts = {};
% % kappats = [];
% % w0 = [0.1867; 0.1761; 0.0623];
% % ksai0 = ksai0_value * rand(nu, 1);
% % 
% % for i_exp = 1:3
% %     p = p_set(i_exp);
% %     windspeed.string1 = p;
% %     windspeed.string2=windspeed.string1;
% %     windspeed.string3=windspeed.string1;
% %     tspan = tspan_set(i_exp);
% % 
% %     out = sim("wind1000_BUS_FOM_ROM");
% %     t = out.tout; 
% %     ksit = squeeze(out.Ksait.Data);
% %     ksit = ksit';
% % 
% %     wt = squeeze(out.Ksait.Data);
% %     wt = wt';
% %     % Set initial condition for next loop
% %     w0 = wt(end, :)';
% %     ksai0 = ksit(end, :)';
% %     % % Dilute data
% %     % t = t(1:sampling2:end-offset2);
% %     % ksit = ksit(1:sampling2:end-offset2, :);
% % 
% %     n_t = size(ksit, 1);
% %     wp = [ksit ones(n_t, 1)*p];
% %     % [wp, m, s] = normalizedata(wp, m, s);
% %     kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;
% % 
% %     ts{i_exp} = t;
% %     kappats = [kappats; kappat];
% % end
% % tss = [ts{1}; ts{2}+40; ts{3}+80];
% % 
% % figure;
% % plot(tss, kappats);
% % hold on
% % % a = xline([0 40 70], ':', {'p = 6.25','p = 9.25', 'p = 7.5'}, 'Color', [17 17 17]/255, 'LineWidth', 1.5); 
% % grid;
% % % ylim([-16, 2]); 
% % % legend("\psi", 'Location', 'southeast');
% % xlabel('t (s)');
% % 
% % load t_444.mat;
% % load yt_444.mat;
% % plot(t,yt);
% % legend("\psi", "y", 'Location', 'southeast');
% 
% % Set wind speed
% p = 9.3;
% % Set sampling settings
% offset2 = 0 * sampling; 
% tspan = 200; % end time of ROM
% sampling2 = 5;
% 
% out = sim("wind1000_BUS_FOM_ROM");
% t = out.tout;
% t = t(1:sampling2:end-offset2);
% ksit = squeeze(out.Ksait.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p];
% % [wp, m, s] = normalizedata(wp, m, s);
% % wp(:, 1) = 0;
% kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;
% 
% plot(Yt_9_3WS.Time(1:50:end-offset2), Yt_9_3WS.Data(1:50:end-offset2));
% % plot(Yt_6_4WS.Time(1:50:end-offset2), Yt_6_4WS.Data(1:50:end-offset2));
% hold on
% plot(t, kappat);
% xlim([80 200]);
% grid on
% legend("FOM", "ROM");

%% Figures 3
power_base = 1e6;
% Set wind speed
p = 6.4;
% Set sampling settings
offset2 = 0 * sampling; 
sampling2 = 50;

t = Ksait_6_4WS.Time(1:sampling2:end-offset2);
ksit = squeeze(Ksait_6_4WS.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p];
kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;

figure;
plot(Yt_6_4WS.Time(1:50:end-offset2), Yt_6_4WS.Data(1:50:end-offset2)./power_base, 'r-', 'LineWidth', 1.2);
hold on
plot(t, kappat./power_base, 'b-.', 'LineWidth', 1.2);
xlim([80 200]);
grid on


% Set wind speed
p = 7.5;
% Set sampling settings
offset2 = 0 * sampling; 
sampling2 = 50;

t = Ksait_7_5WS.Time(1:sampling2:end-offset2);
ksit = squeeze(Ksait_7_5WS.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p];
kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;

% figure;
plot(Yt_7_5WS.Time(1:50:end-offset2), Yt_7_5WS.Data(1:50:end-offset2)./power_base, 'r-', 'LineWidth', 1.2);
hold on
plot(t, kappat./power_base, 'b-.', 'LineWidth', 1.2);
xlim([80 200]);
grid on

% Set wind speed
p = 9.3;
% Set sampling settings
offset2 = 0 * sampling; 
sampling2 = 50;

t = Ksait_9_3WS.Time(1:sampling2:end-offset2);
ksit = squeeze(Ksait_9_3WS.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p];
kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;

% figure;
plot(Yt_9_3WS.Time(1:50:end-offset2), Yt_9_3WS.Data(1:50:end-offset2)./power_base, 'r-', 'LineWidth', 1.2);
hold on
plot(t, kappat./power_base, 'b-.', 'LineWidth', 1.2);
xlim([80 200]);
grid on

% p = 10.3;
% % Set sampling settings
% offset2 = 0 * sampling; 
% sampling2 = 50;
% 
% t = Ksait_10_3WS.Time(1:sampling2:end-offset2);
% ksit = squeeze(Ksait_10_3WS.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p];
% kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;
% 
% % figure;
% plot(Yt_10_3WS.Time(1:50:end-offset2), Yt_10_3WS.Data(1:50:end-offset2)./power_base, 'r-', 'LineWidth', 1.2);
% hold on
% plot(t, kappat./power_base, 'b-.', 'LineWidth', 1.2);
% xlim([80 200]);
% grid on

% p = 5.6;
% % Set sampling settings
% offset2 = 0 * sampling; 
% sampling2 = 50;
% 
% t = Ksait_5_6WS.Time(1:sampling2:end-offset2);
% ksit = squeeze(Ksait_5_6WS.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p];
% kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;
% 
% % figure;
% plot(Yt_5_6WS.Time(1:50:end-offset2), Yt_5_6WS.Data(1:50:end-offset2)./power_base, 'r-', 'LineWidth', 1.2);
% hold on
% plot(t, kappat./power_base, 'b-.', 'LineWidth', 1.2);
% xlim([80 200]);
% grid on


legend("FOM", "ROM");
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 15);
ylabel('Active Power (p.u.)', 'Interpreter', 'latex', 'FontSize', 15);

%% Figures 4
power_base = 1e6;
% Set wind speed
p = 6.4;
% Set sampling settings
offset2 = 0 * sampling; 
sampling2 = 50;

t = Ksait_6_4WS.Time(1:sampling2:end-offset2);
ksit = squeeze(Ksait_6_4WS.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p];
kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;

figure;
% plot(Yt_6_4WS.Time(1:50:end-offset2), (Yt_6_4WS.Data(1:50:end-offset2) - kappat)./power_base, 'r:', 'LineWidth', 1.2);
plot(Yt_6_4WS.Time(1:50:end-offset2), (Yt_6_4WS.Data(1:50:end-offset2) - kappat)./Yt_6_4WS.Data(1:50:end-offset2), 'r:', 'LineWidth', 1.2);
hold on
xlim([80 200]);
grid on
% legend("p = 6.4");


% Set wind speed
p = 7.5;
% Set sampling settings
offset2 = 0 * sampling; 
sampling2 = 50;

t = Ksait_7_5WS.Time(1:sampling2:end-offset2);
ksit = squeeze(Ksait_7_5WS.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p];
kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;

% figure;
% plot(Yt_7_5WS.Time(1:50:end-offset2), (Yt_7_5WS.Data(1:50:end-offset2) - kappat)./power_base, 'g:', 'LineWidth', 1.2);
plot(Yt_7_5WS.Time(1:50:end-offset2), (Yt_7_5WS.Data(1:50:end-offset2) - kappat)./Yt_7_5WS.Data(1:50:end-offset2), 'g:', 'LineWidth', 1.2);
hold on
xlim([80 200]);
grid on
% legend("p = 7.5");


% Set wind speed
p = 9.3;
% Set sampling settings
offset2 = 0 * sampling; 
sampling2 = 50;

t = Ksait_9_3WS.Time(1:sampling2:end-offset2);
ksit = squeeze(Ksait_9_3WS.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p];
kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;

% figure;
% plot(Yt_9_3WS.Time(1:50:end-offset2), (Yt_9_3WS.Data(1:50:end-offset2) - kappat)./power_base, 'b:', 'LineWidth', 1.2);
plot(Yt_9_3WS.Time(1:50:end-offset2), (Yt_9_3WS.Data(1:50:end-offset2) - kappat)./Yt_9_3WS.Data(1:50:end-offset2), 'b:', 'LineWidth', 1.2);
hold on
xlim([80 200]);
grid on
% legend("p = 6.4", "p = 7.5", "p = 9.3");
% xlable('$t$', 'interpreter', 'latex', 'fontsize', 15);
% ylabel('Error (p.u.)', 'Interpreter', 'latex', 'FontSize', 15);

% % Set wind speed
% p = 10.3;
% % Set sampling settings
% offset2 = 0 * sampling; 
% sampling2 = 50;
% 
% t = Ksait_10_3WS.Time(1:sampling2:end-offset2);
% ksit = squeeze(Ksait_10_3WS.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p];
% kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;
% 
% % figure;
% plot(Yt_10_3WS.Time(1:50:end-offset2), (Yt_10_3WS.Data(1:50:end-offset2) - kappat)./power_base, 'm:', 'LineWidth', 1.2);
% hold on
% xlim([80 200]);
% grid on

% % Set wind speed
% p = 5.6;
% % Set sampling settings
% offset2 = 0 * sampling; 
% sampling2 = 50;
% 
% t = Ksait_5_6WS.Time(1:sampling2:end-offset2);
% ksit = squeeze(Ksait_5_6WS.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p];
% kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;
% 
% % figure;
% plot(Yt_5_6WS.Time(1:50:end-offset2), (Yt_5_6WS.Data(1:50:end-offset2) - kappat)./power_base, 'm:', 'LineWidth', 1.2);
% hold on
% xlim([80 200]);
% grid on

legend("p = 6.4", "p = 7.5", "p = 9.3");
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 15);
ylabel('Error (p.u.)', 'Interpreter', 'latex', 'FontSize', 15);