function lambda_w = compute_weights_Milne_SIR(f_tot, x, w, dw, lambda, sigma_J, choice_G, Pmat, Bmat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing weights lambda_w using Milne rule (4th order)    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nw = length(w);

% Define three-points stencil (i+1/4, i+1/2, i+3/4)

w1 = 1/4 * (3*w(1:Nw-1) + w(2:Nw));   % First point
w2 = 1/2 * (w(1:Nw-1) + w(2:Nw));     % Middle point
w3 = 1/4 * (w(1:Nw-1) + 3*w(2:Nw));   % Third point 


% Define weights D_w

D1 = diffusion_SIR(f_tot, x, w1, w, sigma_J, Bmat);
D2 = diffusion_SIR(f_tot, x, w2, w, sigma_J, Bmat);
D3 = diffusion_SIR(f_tot, x, w3, w, sigma_J, Bmat);


% Define partial derivatives of diffusion term

dD1 = partial_diffusion_SIR(f_tot, x, w1, w, sigma_J, Bmat);
dD2 = partial_diffusion_SIR(f_tot, x, w2, w, sigma_J, Bmat);
dD3 = partial_diffusion_SIR(f_tot, x, w3, w, sigma_J, Bmat);


% Define drift terms

K1 = drift_SIR(f_tot, x, w1, w, lambda, choice_G, Pmat, Bmat);
K2 = drift_SIR(f_tot, x, w2, w, lambda, choice_G, Pmat, Bmat);
K3 = drift_SIR(f_tot, x, w3, w, lambda, choice_G, Pmat, Bmat);


% Define weights C_w

C1 = K1 + dD1;
C2 = K2 + dD2;
C3 = K3 + dD3;


% Define weights lambda_w using Milne rule

lambda_w = dw/3 * ( 2*C1./D1 - C2./D2 + 2*C3./D3 );    


end