function [rho, mean] = compute_moments_popularity(f, v, dv)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Function computing density and mean popularity of distribution  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute macroscopic density

rho = sum(f) * dv;


% Compute mean popularity

mean = 1/rho * sum(f .* v) * dv;


end