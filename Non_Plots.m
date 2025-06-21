%% Figures
power_base = 1e6;
% Set wind speed
p = 7.5;
% Set sampling settings
offset2 = 0 * sampling; 
sampling2 = 50;

Ksait = out.Ksait;
Yt = out.Yt;

t = Ksait.Time(1:sampling2:end-offset2);
ksit = squeeze(Ksait.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit p.*ones(n_t, 1)];
kappat = RBFbasisnD(wp, centers, sigmas) * coeffs;

% Create plots.
tt = tiledlayout(2,1);
ax1 = nexttile;
% figure;
plot(ax1, Yt.Time(1:50:end-offset2), Yt.Data(1:50:end-offset2)./power_base, 'r-', 'LineWidth', 1.2);
hold on
plot(ax1, t, kappat./power_base, 'b-.', 'LineWidth', 1.2);
xlim([82 142]);
grid on
legend("FOM", "ROM");
ylabel('Active Power (p.u.)', 'Interpreter', 'latex', 'FontSize', 15);

ax2 = nexttile;
plot(ax2, Yt.Time(1:50:end-offset2), abs(Yt.Data(1:50:end-offset2) - kappat)./abs(Yt.Data(1:50:end-offset2)), 'g:', 'LineWidth', 1.2);
hold on
xlim([82 142]);
grid on
ylabel('Error (p.u.)', 'Interpreter', 'latex', 'FontSize', 15);

% Link the axes
linkaxes([ax1,ax2],'x');
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 15);

% Move plots closer together
xticklabels(ax1,{});
tt.TileSpacing = 'compact';