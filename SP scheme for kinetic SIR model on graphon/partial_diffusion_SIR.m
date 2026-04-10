function d_D = partial_diffusion_SIR(f_tot, x, w_eval, w, sigma_J, Bmat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing partial derivative of diffusion term    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dy = x(2) - x(1);
dv = w(2) - w(1);

rho_x = sum(f_tot, 2) * dv;   % Integral in v of f_tot          

By_f = Bmat * rho_x * dy;   % Weighted integral in y depending on pre-computed graphon


% Compute partial derivative of diffusion term on shifted grid w_eval

d_D = -sigma_J * (By_f * w_eval(:).');


end