function plot_class_snapshots(w, fJ_over_time, ix, snap, t_snap, BJ, Jlabel, STYLE)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function plotting selected snapshots of f_J, optionally including steady state    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Snapshot colors

C_ocra        = [0.90 0.62 0.00];
C_green       = [0.00 0.55 0.25];
C_purple      = [0.50 0.20 0.70];
C_blue        = [0.00 0.45 0.74];
C_orange      = [0.85 0.33 0.10];
C_siena       = [0.60 0.24 0.10];
C_orange_soft = [0.92 0.60 0.40];
C_t0          = [0.35 0.40 0.45];


% Steady states color

C_red         = [0.85 0.20 0.20];


% Colors for fixed times

Tmap = [
    10   C_green
    20   C_siena
    30   C_purple
    35   C_ocra
    110  C_blue
];

fig = new_figure(STYLE);
ax = gca;
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'off');

set(ax, 'FontSize', STYLE.fsAxes, ...
    'LineWidth', STYLE.lwAxes, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'off', 'YMinorTick', 'off');


% Plot snapshots

for s = 1:numel(snap)

    k = snap(s);
    y = squeeze(fJ_over_time(ix, :, k));

    t_req = t_snap(s);

    if abs(t_req) < 1e-12

        c = C_t0;

    else

        row = find(abs(Tmap(:, 1) - t_req) < 1e-12, 1);

        if isempty(row)
            
            palette = [C_orange; C_green; C_purple; C_blue; C_orange_soft];

            c = palette(1 + mod(s-1, size(palette, 1)), :);

        else

            c = Tmap(row, 2:4);

        end

    end

    plot(ax, w, y, 'LineWidth', STYLE.lw, 'Color', c, 'LineStyle','-');

end


% Steady state flag, true if available

hasSteady = (nargin >= 7) && ~isempty(BJ);

if hasSteady

    plot(ax, w, BJ(ix, :), 'LineWidth', STYLE.lwSteady+1, ...
         'Color', C_red, 'LineStyle', '--');
end

xlabel(ax, '$w$', 'FontSize', STYLE.fsLabel, 'Interpreter', 'latex');
ylabel(ax, sprintf('$f_{%s}(x^\\star,w)$', Jlabel), ...
       'FontSize', STYLE.fsLabel, 'Interpreter', 'latex');

xlim(ax, [-1 1]);
xticks(ax, -1:0.5:1);


% Legend

if hasSteady

    leg = cell(1, numel(snap)+1);

else

    leg = cell(1,numel(snap));

end

for s = 1:numel(snap)

    leg{s} = sprintf('$t=%.2f$', t_snap(s));

end

if hasSteady

    leg{end} = 'Steady';

end

lg = legend(ax, leg, 'Location', 'northwest', 'Interpreter', 'latex');
set(lg, 'Box', 'on', 'FontSize', STYLE.fsLegend);

if isprop(lg,'ItemTokenSize')

    lg.ItemTokenSize = [65 25];

end

hold(ax, 'off');


end