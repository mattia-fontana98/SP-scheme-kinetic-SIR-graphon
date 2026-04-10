function h_steady = compute_steady_state_popularity(v, dv, mu, zeta, theta, F_inf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function computing steady state (inverse Gamma distribution)    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inverse Gamma parameters

alpha = 1 + 2 * mu/zeta;
beta  = 2 * theta * F_inf/zeta;


% Avoid singularity at v = 0

v_safe = v;
v_safe(1) = dv/2;


% Compute steady state

h_steady = (beta^alpha / gamma(alpha)) .* (v_safe.^(-alpha-1)) .* exp(-beta ./ v_safe);


end