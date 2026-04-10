function lambda_v = compute_weights_NC6_popularity(F, v, mu, theta, zeta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing weights lambda_v using Newton-Cotes rule (6th order)    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nv = length(v);
dv = v(2) - v(1);


% Define five-point stencil (i+1/6, i+1/3, i+1/2, i+2/3, i+5/6)

v1 = 1/6 * (v(2:Nv) + 5*v(1:Nv-1));   % First point
v2 = 1/3 * (v(2:Nv) + 2*v(1:Nv-1));   % Second point
v3 = 1/2 * (v(1:Nv-1) + v(2:Nv));     % Middle point 
v4 = 1/3 * (2*v(2:Nv) + v(1:Nv-1));   % Fourth point
v5 = 1/6 * (5*v(2:Nv) + v(1:Nv-1));   % Fifth point 


% Define weights D_v

D1 = diffusion_popularity(v1, zeta);   
D2 = diffusion_popularity(v2, zeta);   
D3 = diffusion_popularity(v3, zeta);
D4 = diffusion_popularity(v4, zeta);   
D5 = diffusion_popularity(v5, zeta);


% Define partial derivatives of diffusion term

dD1 = partial_diffusion_popularity(v1, zeta);
dD2 = partial_diffusion_popularity(v2, zeta);
dD3 = partial_diffusion_popularity(v3, zeta);
dD4 = partial_diffusion_popularity(v4, zeta);
dD5 = partial_diffusion_popularity(v5, zeta);


% Define drift terms

K1 = drift_popularity(F, v1, mu, theta);
K2 = drift_popularity(F, v2, mu, theta);
K3 = drift_popularity(F, v3, mu, theta);
K4 = drift_popularity(F, v4, mu, theta);
K5 = drift_popularity(F, v5, mu, theta);


% Define weights C_v

C1 = K1 + dD1;
C2 = K2 + dD2;
C3 = K3 + dD3;
C4 = K4 + dD4;
C5 = K5 + dD5;


% Define weights lambda_v using Newton-Cotes rule with 5 points

lambda_v = dv/20 * ( 11*C1./D1 - 14*C2./D2 + 26*C3./D3 - 14*C4./D4 + 11*C5./D5 );


end