
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Plot script for kinetic SIR model with graphon    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear

clc

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 18);



%% Load data


load(fullfile('..', 'Output_kineticSIR.mat'));



%% Plot style


STYLE = setup_plot_style();



%% User switches


% Control plot (0 --> Deactivate plot / 1 --> Activate plot)

plot_final_surfaces = 1;


% Plot time evolution of rho_J

plot_compartment_masses = 1;


% Plot time evolution of effective reproduction number 

plot_reproduction_number = 1;


% Plot of errors, when available

plot_errors = 1;


% Plot compartmental snapshots

plot_snapshots = 1;



%% User parameters for snapshots


ix_S = 3;                % x-index used in snapshot plot of susceptibles
ix_I = 11;               % x-index used in snapshot plot of infected
ix_R = 19;               % x-index used in snapshot plot of recovered

iw    = 2:length(w)-1;   % Use 1:length(w) to include boundaries


% Snapshot times

t_snap_S = [0 10 110];
t_snap_I = [10 30 35 110];
t_snap_R = [10 20 30 110];



%% Final surfaces


if plot_final_surfaces == 1

    new_figure(STYLE);
    surf(w, x, fS);
    set(gca, 'FontSize', STYLE.fsAxes, 'LineWidth', STYLE.lwAxes, 'TickLabelInterpreter', 'latex');
    xlabel('$w$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    ylabel('$x$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    zlabel('$f_S(T,x,w)$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);

    new_figure(STYLE);
    surf(w, x, fS_over_time(:, :, 1));
    set(gca, 'FontSize', STYLE.fsAxes, 'LineWidth', STYLE.lwAxes, 'TickLabelInterpreter', 'latex');
    xlabel('$w$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    ylabel('$x$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    zlabel('$f_S^{in}(x,w)$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);

    new_figure(STYLE);
    surf(w, x, fR);
    set(gca, 'FontSize', STYLE.fsAxes, 'LineWidth', STYLE.lwAxes, 'TickLabelInterpreter', 'latex');
    xlabel('$w$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    ylabel('$x$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    zlabel('$f_R(T,x,w)$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);

    new_figure(STYLE);
    surf(w, x, fI);
    set(gca, 'FontSize', STYLE.fsAxes, 'LineWidth', STYLE.lwAxes, 'TickLabelInterpreter', 'latex');
    xlabel('$w$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    ylabel('$x$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    zlabel('$f_{I}(T,x,w)$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);

end



%% Masses of compartments over time


if plot_compartment_masses == 1

    colS = [0.00 0.45 0.70];
    colI = [0.90 0.62 0.00];
    colR = [0.00 0.62 0.45];

    new_figure(STYLE);
    plot(time, rhoS_over_time, 'LineWidth', 5, 'Color', colS); hold on;
    plot(time, rhoI_over_time, 'LineWidth', 5, 'Color', colI);
    plot(time, rhoR_over_time, 'LineWidth', 5, 'Color', colR);
    hold off

    set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLabelInterpreter', 'latex');
    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off');

    xlabel('$t$', 'FontSize', STYLE.fsLabel);
    ylabel('$\rho_J(t)$', 'FontSize', STYLE.fsLabel);

    lg = legend({'$S$','$I$','$R$'}, 'Location', 'southwest', 'Interpreter', 'latex');
    set(lg, 'Box', 'on', 'FontSize', 35);
    lg.ItemTokenSize = [55 18];

end



%% Effective reproduction number


if plot_reproduction_number == 1 && alpha == 1

    new_figure(STYLE);
    plot(time, R_over_time, 'LineWidth', 3);
    hold on
    yline(1, '--');
    hold off

    set(gca, 'FontSize', STYLE.fsAxes, 'LineWidth', STYLE.lwAxes, 'TickLabelInterpreter', 'latex');
    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);
    ylabel('$R_{\mathrm{eff}}(t)$', 'Interpreter', 'latex', 'FontSize', STYLE.fsLabel);

end



%% Errors


if plot_errors == 1 && exist('error_S', 'var') && exist('error_I', 'var') && exist('error_R', 'var')

    has_error_S = any(isfinite(error_S) & (error_S > 0));
    has_error_I = any(isfinite(error_I) & (error_I > 0));
    has_error_R = any(isfinite(error_R) & (error_R > 0));

    if has_error_S 

        colS = [0.00 0.45 0.70];
        colI = [0.90 0.62 0.00];
        colR = [0.00 0.62 0.45];

        new_figure(STYLE);

        leg = {};

            semilogy(time, error_S, 'LineWidth', 5, 'Color', colS); hold on
            leg{end+1} = '$S$';
            semilogy(time, error_I, 'LineWidth', 5, 'Color', colI);
            leg{end+1} = '$I$';

            semilogy(time, error_R, 'LineWidth', 5, 'Color', colR);
            leg{end+1} = '$R$';
        hold off

        set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLabelInterpreter', 'latex');
        set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off');

        xlabel('$t$', 'FontSize', STYLE.fsLabel);
        ylabel('$\|f_J(t,\cdot)-f_J^\infty\|_{L^1(\Omega\times\mathcal{I})}$', 'FontSize', STYLE.fsLabel);

        lg = legend(leg, 'Location', 'southwest', 'Interpreter', 'latex');
        set(lg, 'Box', 'on', 'FontSize', 35);
        lg.ItemTokenSize = [55 18];

    end

end



%% Snapshots in w at fixed x-index


if plot_snapshots == 1

    % Map requested times to nearest available time indices

    snap_S = get_snap_indices(time, t_snap_S);
    snap_I = get_snap_indices(time, t_snap_I);
    snap_R = get_snap_indices(time, t_snap_R);


    % Steady states, if available

    if exist('steady_state_S', 'var')

        BJ_S = steady_state_S(:, iw);

    else

        BJ_S = [];

    end

    if exist('steady_state_I', 'var')

        BJ_I = steady_state_I(:, iw);

    else

        BJ_I = [];

    end

    if exist('steady_state_R', 'var')

        BJ_R = steady_state_R(:, iw);

    else

        BJ_R = [];

    end

plot_class_snapshots(w(iw), fS_over_time(:, iw, :), 3, snap_S, t_snap_S, BJ_S, 'S', STYLE);
plot_class_snapshots(w(iw), fI_over_time(:, iw, :), 11, snap_I, t_snap_I, BJ_I, 'I', STYLE);
plot_class_snapshots(w(iw), fR_over_time(:, iw, :), 19, snap_R, t_snap_R, BJ_R, 'R', STYLE);

end
