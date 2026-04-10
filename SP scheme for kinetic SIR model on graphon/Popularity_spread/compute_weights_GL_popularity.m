function lambda_v = compute_weights_GL_popularity(F, v, mu, theta, zeta, m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Function computing weights lambda_v using Gauss–Legendre quadrature of order m    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nv = length(v);

lambda_v = zeros(1, Nv-1);

for i = 1:(Nv-1)

    % Gauss-Legendre nodes and weights on [v_i, v_{i+1}]

    [xj, wj] = lgwt(m, v(i), v(i+1));
    xj = xj(:).';
    wj = wj(:).';


    % Evaluate drift and diffusion terms at quadrature nodes

    Ck = drift_popularity(F, xj, mu, theta) + partial_diffusion_popularity(xj, zeta);

    Dk = diffusion_popularity(xj, zeta);


    % Quadrature approximation of ∫_[0,L] C(v)/D(v) dv

    lambda_v(i) = sum( wj .* (Ck ./ Dk) );
    
end


end