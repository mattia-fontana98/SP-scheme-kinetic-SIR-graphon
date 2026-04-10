function [steady_state, rho_P, m_P] = compute_steady_state_SIR(f_tot, f_in, x, dx, w, dw, lambda, sigma_J, Px)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing steady state (Beta distribution)    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[W, ~] = meshgrid(w, x);

[Nx, Nw] = size(W);

nu = sigma_J / lambda;

steady_state = zeros(Nx, Nw);


% Compute weighted mass and mean

rho_P = sum(Px .* f_tot, 'all') * dw * dx;

m_P = sum(W .* Px .* f_tot, 'all') * dw * dx;

cutoff1 = 1e-10;

cutoff2 = 1e-10;


% Compute steady state

for i = 1:Nx

    fw = @(w) (1 + w).^(-1 + Px(i) / nu * (rho_P + m_P)) .* (1 - w).^(-1 + Px(i) / nu * (rho_P - m_P));

    C_norm = (sum(f_in(i, :)) * dw) / integral(fw, -1+cutoff1, 1-cutoff1);

    steady_state(i, :) = C_norm * (1 + w + cutoff2 * (w == -1)).^(-1 + Px(i) / nu * (rho_P + m_P)) .* (1 - w + cutoff2 * (w == 1)).^(-1 + Px(i) / nu * (rho_P - m_P));

end


end