function D = diffusion_SIR(f_tot, x, w_eval, w, sigma_J, Bmat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing diffusion term    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dy = x(2) - x(1);
dv = w(2) - w(1);

rho_x = sum(f_tot, 2) * dv;   % Integral in v of f_tot            

By_f = Bmat * rho_x * dy;   % Weighted integral in y depending on pre-computed graphon            


% Compute diffusion term on shifted grid w_eval

D = (sigma_J/2) * (By_f * (1 - w_eval(:).'.^2));


end