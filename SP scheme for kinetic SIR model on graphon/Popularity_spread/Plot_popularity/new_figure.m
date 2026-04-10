function fig = new_figure(STYLE)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function creating new figure in fixed dimensions and position    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fig = figure('Color', 'w');

set(fig,'Units', 'centimeters');

set(fig, 'Position', [2 2 STYLE.Wcm STYLE.Hcm]);

set(fig, 'PaperUnits', 'centimeters');

set(fig, 'PaperPositionMode', 'manual');

set(fig, 'PaperPosition', [0 0 STYLE.Wcm STYLE.Hcm]);

set(fig, 'PaperSize', [STYLE.Wcm STYLE.Hcm]);


end