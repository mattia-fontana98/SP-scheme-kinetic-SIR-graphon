function [delta_w, C_w, D_w] = compute_weights_SIR(f_tot, x, w, dw, lambda, sigma_J, choice_G, Pmat, Bmat, scheme_w)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing weights delta_w, C_w and D_w for Chang-Cooper scheme    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nw = length(w);


% Shift to intermediate nodes (i+1/2,j)

w_mid = 1/2 * (w(1:Nw-1) + w(2:Nw));


% Compute weights D_w

D_w = diffusion_SIR(f_tot, x, w_mid, w, sigma_J, Bmat);


% Compute weights C_w

C_w = drift_SIR(f_tot, x, w_mid, w, lambda, choice_G, Pmat, Bmat) + partial_diffusion_SIR(f_tot, x, w_mid, w, sigma_J, Bmat);


% Compute weights lambda_w

lambda_midpoint = dw * C_w ./ D_w;

if scheme_w == 1   % Midpoint rule (2nd order)

    lambda_w = dw * C_w ./ D_w;

elseif scheme_w == 2   % Milne rule (4th order)

    lambda_w = compute_weights_Milne_SIR(f_tot, x, w, dw, lambda, sigma_J, choice_G, Pmat, Bmat);

elseif scheme_w == 3   % Newton-Cotes rule (6th order)

    lambda_w = compute_weights_NC6_SIR(f_tot, x, w, dw, lambda, sigma_J, choice_G, Pmat, Bmat);

elseif scheme_w == 4   % Gauss-Legendre rule (8th order)

    lambda_w = compute_weights_GL_SIR(f_tot, x, w, lambda, sigma_J, choice_G, Pmat, Bmat, 8);
end


% Compute weights delta_w

tol = 1e-10;

cutoff_lambda = (abs(lambda_w) < tol);

delta_w = ( 1 ./ lambda_midpoint + 1 ./ (1 - exp(lambda_w)) ) .* (1 - cutoff_lambda) + cutoff_lambda / 2;


end