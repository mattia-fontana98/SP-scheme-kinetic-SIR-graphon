
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Plot script for popularity spread   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear
clc



%% Load output data


load(fullfile('..', 'Output_popularity.mat'));



%% Plot style


STYLE = setup_plot_style();



%% Variables used in plots

tH = time(:);
h_over_time = h_save;



%% User switches


plot_snapshots = 1;
plot_errors    = 1;



%% User parameters


t_req = [0 0.1 100];



%% Plot snapshots of h(t,v)


if plot_snapshots == 1

    plot_h_snapshots(v, tH, h_over_time, h_steady, t_req, STYLE, ...
                     STYLE.snapColors, STYLE.snapStyles, STYLE.steadyColor, STYLE.steadyStyle);

end



%% Plot convergence errors


if plot_errors == 1

    fig = new_figure(STYLE);
    ax = axes('Parent', fig);
    
    yyaxis(ax, 'left')
    p1 = semilogy(ax, tH, max(errL1,1e-16), 'LineWidth', STYLE.lw);
    
    yyaxis(ax, 'right')
    p2 = semilogy(ax, tH, max(eF,1e-16), 'LineWidth', STYLE.lw);
    
    set(ax, 'FontSize', STYLE.fsAxes, 'LineWidth', STYLE.lwAxes, 'TickLabelInterpreter', ...
        'latex', 'XMinorTick', 'off', 'YMinorTick', 'off');

    ax.YAxis(1).MinorTickValues = [];
    ax.YAxis(2).MinorTickValues = [];

    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);

    lg = legend([p1 p2], {'$\|h-h^\infty\|_{L^1(R_+)}$','$|\mathcal{F}-\mathcal{F}^\infty|$'}, ...
                'Interpreter', 'latex', 'Location', 'northeast');

    set(lg, 'Box', 'on', 'FontSize', STYLE.fsLegend);

    ax.YAxis(1).Color = p1.Color;
    ax.YAxis(2).Color = p2.Color;


end
