function STYLE = setup_plot_style()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function defining style of plots    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Figure size

STYLE.Wcm = 30;
STYLE.Hcm = 20;


% Font sizes

STYLE.fsAxes   = 30;
STYLE.fsLabel  = 35;
STYLE.fsTitle  = 25;
STYLE.fsLegend = 30;


% Line widths

STYLE.lw       = 5;
STYLE.lwSteady = 5.5;
STYLE.lwAxes   = 2;


% Class colors

STYLE.colS = [0.00 0.45 0.74];
STYLE.colI = [0.85 0.20 0.20];
STYLE.colR = [0.00 0.40 0.20];


% Extra colors

STYLE.C_green  = [0.00 0.55 0.25];
STYLE.C_purple = [0.50 0.20 0.70];
STYLE.C_blue   = [0.00 0.45 0.74];
STYLE.C_red    = [0.85 0.20 0.20];


% Snapshot styling

STYLE.snapColors = [STYLE.C_green; STYLE.C_purple; STYLE.C_blue];
STYLE.snapStyles = {'-', '-', '-'};


% Steady state styling

STYLE.steadyColor = STYLE.C_red;
STYLE.steadyStyle = '--';


% Global LaTeX settings

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');


end