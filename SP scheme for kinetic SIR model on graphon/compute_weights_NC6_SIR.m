function lambda_w = compute_weights_NC6_SIR(f_tot, x, w, dw, lambda, sigma_J, choice_G, Pmat, Bmat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing weights lambda_w using Newton-Cotes rule (6th order)    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nw = length(w);


% Define five-points stencil (i+1/6, i+1/3, i+1/2, i+2/3, i+5/6)

w1 = 1/6 * (w(2:Nw) + 5*w(1:Nw-1));   % First point
w2 = 1/3 * (w(2:Nw) + 2*w(1:Nw-1));   % Second point
w3 = 1/2 * (w(1:Nw-1) + w(2:Nw));     % Middle point 
w4 = 1/3 * (2*w(2:Nw) + w(1:Nw-1));   % Fourth point
w5 = 1/6 * (5*w(2:Nw) + w(1:Nw-1));   % Fifth point 


% Define weights D_w

D1 = diffusion_SIR(f_tot, x, w1, w, sigma_J, Bmat);
D2 = diffusion_SIR(f_tot, x, w2, w, sigma_J, Bmat);
D3 = diffusion_SIR(f_tot, x, w3, w, sigma_J, Bmat);
D4 = diffusion_SIR(f_tot, x, w4, w, sigma_J, Bmat);
D5 = diffusion_SIR(f_tot, x, w5, w, sigma_J, Bmat);


% Define partial derivatives of diffusion term

dD1 = partial_diffusion_SIR(f_tot, x, w1, w, sigma_J, Bmat);
dD2 = partial_diffusion_SIR(f_tot, x, w2, w, sigma_J, Bmat);
dD3 = partial_diffusion_SIR(f_tot, x, w3, w, sigma_J, Bmat);
dD4 = partial_diffusion_SIR(f_tot, x, w4, w, sigma_J, Bmat);
dD5 = partial_diffusion_SIR(f_tot, x, w5, w, sigma_J, Bmat);


% Define drift terms

K1 = drift_SIR(f_tot, x, w1, w, lambda, choice_G, Pmat, Bmat);
K2 = drift_SIR(f_tot, x, w2, w, lambda, choice_G, Pmat, Bmat);
K3 = drift_SIR(f_tot, x, w3, w, lambda, choice_G, Pmat, Bmat);
K4 = drift_SIR(f_tot, x, w4, w, lambda, choice_G, Pmat, Bmat);
K5 = drift_SIR(f_tot, x, w5, w, lambda, choice_G, Pmat, Bmat);


% Define weights C_w

C1 = K1 + dD1;
C2 = K2 + dD2;
C3 = K3 + dD3;
C4 = K4 + dD4;
C5 = K5 + dD5;


% Define weights lambda_w using Newton-Cotes rule with 5 points

lambda_w = dw/20 * ( 11*C1./D1 - 14*C2./D2 + 26*C3./D3 - 14*C4./D4 + 11*C5./D5 );


end
