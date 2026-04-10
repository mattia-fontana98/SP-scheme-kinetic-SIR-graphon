function [delta_v, C_v, D_v] = compute_weights_popularity(F, v, dv, mu, zeta, theta, scheme_v)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing weights delta_v, C_v and D_v for Chang-Cooper scheme    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nv = length(v);


% Shift domain to intermediate nodes i+1/2
    
v_mid = ( v(1:Nv-1) + v(2:Nv)) / 2;
    
    
% Compute weights D_v

D_v = diffusion_popularity(v_mid, zeta);


% Compute weights C_v

C_v = drift_popularity(F, v_mid, mu, theta) + partial_diffusion_popularity(v_mid, zeta);

lambda_midpoint = dv * C_v ./ D_v; 


% Compute weights lambda_v

if scheme_v == 1   % Midpoint (2nd order)
        
    lambda_v = dv * C_v ./ D_v; 

elseif scheme_v == 2   % Milne (4th order)

    lambda_v = compute_weights_Milne_popularity(F, v, mu, theta, zeta);

elseif scheme_v == 3   % Newton-Cotes (6th order)

    lambda_v = compute_weights_NC6_popularity(F, v, mu, theta, zeta);

elseif scheme_v == 4   % Gauss-Legendre (8th order)

    lambda_v = compute_weights_GL_popularity(F, v, mu, theta, zeta, 8);

end


% Compute weights delta_v

tol = 1e-10;

cutoff_lambda = (abs(lambda_v) < tol);

delta_v = ( 1 ./ lambda_midpoint + 1 ./ (1 - exp(lambda_v)) ) .* (1 - cutoff_lambda) + cutoff_lambda / 2;


end