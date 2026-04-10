function plot_h_snapshots(v, tH, h_over_time, h_steady, t_req, STYLE, snapColors, snapStyles, steadyColor, steadyStyle)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function plotting selected snapshots of h, together with steady state    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


new_figure(STYLE);
ax = gca;
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'off');

set(ax, 'FontSize', STYLE.fsAxes, 'LineWidth', STYLE.lwAxes, 'TickLabelInterpreter', ...
    'latex', 'XMinorTick', 'off', 'YMinorTick', 'off');


% Map requested times to nearest available time indices

ns  = numel(t_req);
idx = zeros(ns, 1);

for s = 1:ns

    [~, idx(s)] = min(abs(tH - t_req(s)));

end


% Plot snapshots

for s = 1:numel(idx)

    k = idx(s);
    h = h_over_time(:, k);

    c  = snapColors(1 + mod(s-1,size(snapColors, 1)), :);
    ls = snapStyles{1 + mod(s-1,numel(snapStyles))};

    plot(ax, v, h, 'LineWidth', STYLE.lw, 'Color', c, 'LineStyle', ls);

end


% Plot steady state

plot(ax, v, h_steady, 'LineWidth', STYLE.lwSteady, 'Color', steadyColor, ...
     'LineStyle', steadyStyle);

xlabel(ax, '$v$', 'FontSize', STYLE.fsLabel, 'Interpreter', 'latex');
ylabel(ax, '$h(t,v)$', 'FontSize', STYLE.fsLabel, 'Interpreter', 'latex');

xlim(ax,[v(1) v(end)]);


% Legend

leg = cell(1, numel(idx)+1);

for s = 1:numel(idx)

    leg{s} = sprintf('$t=%.3f$', tH(idx(s)));

end

leg{end} = '$h_\infty$';

lg = legend(ax, leg, 'Location', 'northeast', 'Interpreter', 'latex');
set(lg, 'Box', 'on', 'FontSize', STYLE.fsLegend);

if isprop(lg,'ItemTokenSize')

    lg.ItemTokenSize = [65 25];

end

hold(ax, 'off');


end