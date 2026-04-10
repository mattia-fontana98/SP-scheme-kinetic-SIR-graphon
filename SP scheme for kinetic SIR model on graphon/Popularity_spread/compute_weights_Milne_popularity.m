function lambda_v = compute_weights_Milne_popularity(F, v, mu, theta, zeta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing weights lambda_v using Milne rule (4th order)    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nv = length(v);
dv = v(2) - v(1);


% Define three-point stencil (i+1/4, i+1/2, i+3/4)

v1 = 1/4 * (3*v(1:Nv-1) + v(2:Nv));   % First point
v2 = 1/2 * (v(1:Nv-1) + v(2:Nv));     % Middle point
v3 = 1/4 * (v(1:Nv-1) + 3*v(2:Nv));   % Third point 


% Define weights D_v

D1 = diffusion_popularity(v1, zeta);
D2 = diffusion_popularity(v2, zeta);
D3 = diffusion_popularity(v3, zeta);


% Define partial derivatives of diffusion term

dD1 = partial_diffusion_popularity(v1, zeta);
dD2 = partial_diffusion_popularity(v2, zeta);
dD3 = partial_diffusion_popularity(v3, zeta);


% Define drift terms

K1 = drift_popularity(F, v1, mu, theta);
K2 = drift_popularity(F, v2, mu, theta);
K3 = drift_popularity(F, v3, mu, theta);


% Define weights C_v

C1 = K1 + dD1;
C2 = K2 + dD2;
C3 = K3 + dD3;


% Define weights lambda_v using Milne rule

lambda_v = dv/3 * ( 2*C1./D1 - C2./D2 + 2*C3./D3 );    


end