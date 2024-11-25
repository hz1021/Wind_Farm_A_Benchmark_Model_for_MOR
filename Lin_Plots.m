%% Bode plot
xlim_min = 0; xlim_max = 6.1;
xlim_bode = logspace(xlim_min, xlim_max, 1000);
ylim_min = -10; ylim_max= 10;

interpolated = [imag(eigenS), imag(eigenQ)];
interpolated_filter = interpolated(interpolated ~= 0);

[mag_FOM,phase_FOM, wout_FOM] = bode(sys, xlim_bode);
[mag_ROM,phase_ROM, wout_ROM] = bode(sysr, xlim_bode);
[mag_int,phase_int, wout_int] = bode(sys, interpolated_filter);

% wrap the phase of ROM
phase_ROM =  phase_ROM(:);
min_phase = min(phase_FOM(:)) - 90 - 5.5;
phase_ROM = phase_ROM - min_phase;
phase_ROM  = mod(phase_ROM, 360);
phase_ROM = phase_ROM + min_phase;

mag_FOM = 20 * log10(mag_FOM(:)); 
mag_ROM = 20 * log10(mag_ROM(:));
mag_int = 20 * log10(mag_int(:));

% Create plots.
t = tiledlayout(2,1);
ax1 = nexttile;

semilogx(ax1, wout_FOM, squeeze(mag_FOM),'LineWidth', 1.3, 'Color', 'b');
hold on;
semilogx(ax1, wout_ROM, squeeze(mag_ROM),'LineWidth', 1.1, 'Color', 'r', 'LineStyle', '--');
hold on;
semilogx(ax1, interpolated_filter, squeeze(mag_int), 'o', 'MarkerSize',5, 'MarkerEdgeColor','black', 'MarkerFaceColor','Green');
xlim([10^xlim_min 10^xlim_max]);
% ylim([ylim_min ylim_max]);
% grid on
xline(linspace(1,10,10), ':');
xline(linspace(10,100,10), ':');
xline(linspace(100,1000,10), ':');
xline(linspace(1e3,1e4,10), ':');
xline(linspace(1e4,1e5,10), ':');
xline(linspace(1e5,1e6,10), ':');
yline([110 120 130 140 150 160 170], ':');
ylabel('Magnitude (dB)', 'Interpreter', 'latex', 'FontSize', 15);
% xlabel('Frequency $\omega$', 'Interpreter', 'latex', 'FontSize', 15);
% title('Bode Diagram', 'Interpreter', 'latex', 'FontSize', 15);
legend('FOM', 'ROM', 'Interpreter', 'latex', 'FontSize', 12);

ax2 = nexttile;
semilogx(ax2, wout_FOM, squeeze(phase_FOM),'LineWidth', 1.3, 'Color', 'b');
hold on;
semilogx(ax2, wout_ROM, squeeze(phase_ROM),'LineWidth', 1.1, 'Color', 'r', 'LineStyle', '--');
hold on;
semilogx(ax2, interpolated_filter,squeeze(phase_int), 'o', 'MarkerSize',5, 'MarkerEdgeColor','black', 'MarkerFaceColor','Green');
xlim([10^xlim_min 10^xlim_max]);
xline(linspace(1,10,10), ':');
xline(linspace(10,100,10), ':');
xline(linspace(100,1000,10), ':');
xline(linspace(1e3,1e4,10), ':');
xline(linspace(1e4,1e5,10), ':');
xline(linspace(1e5,1e6,10), ':');
yline([100 150 200 250], ':');
ylabel('Phase (deg)', 'Interpreter', 'latex', 'FontSize', 15);
% xlabel('Frequency $\omega$', 'Interpreter', 'latex', 'FontSize', 15);
% title('Bode Diagram', 'Interpreter', 'latex', 'FontSize', 15);
legend('FOM', 'ROM', 'Interpreter', 'latex', 'FontSize', 12);

% Link the axes
linkaxes([ax1,ax2],'x');

% Add shared title and axis labels
title(t, 'Bode Diagram', 'Interpreter', 'latex', 'FontSize', 15);
xlabel(t, 'Frequency $\omega$', 'Interpreter', 'latex', 'FontSize', 15);

% Move plots closer together
xticklabels(ax1,{});
t.TileSpacing = 'compact';