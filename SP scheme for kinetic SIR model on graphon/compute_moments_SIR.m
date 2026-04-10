function [rho, mean] = compute_moments_SIR(f, x, dx, w, dw)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing density and mean opinion of distributions   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[W, ~] = meshgrid(w, x);


% Compute macroscopic density

rho = sum(f, 'all') * dw * dx;


% Compute mean opinion

mean = 1/(max(rho, 1e-14)) * sum(W .* f, 'all') * dw * dx;


end