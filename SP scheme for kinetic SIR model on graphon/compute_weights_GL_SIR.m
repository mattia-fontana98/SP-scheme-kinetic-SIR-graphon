function lambda_w = compute_weights_GL_SIR(f_tot, x, w, lambda, sigma_J, choice_G, Pmat, Bmat, m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Function computing weights lambda_w using Gauss–Legendre quadrature of order m    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Nw = numel(w);
    Nx = numel(x);

    lambda_w = zeros(Nx, Nw-1);

    for i = 1:(Nw-1)

        % Gauss-Legendre nodes and weights on [w_i, w_{i+1}]

        [wj_nodes, wj_weights] = lgwt(m, w(i), w(i+1));
        wj_nodes   = wj_nodes(:).';    
        wj_weights = wj_weights(:).';  


        acc = zeros(Nx, 1);

        for k = 1:m
            
            wk = wj_nodes(k);


            % Evaluate drift and diffusion terms at the quadrature nodes

            Dk = diffusion_SIR(f_tot, x, wk, w, sigma_J, Bmat); 

            Ck = drift_SIR(f_tot, x, wk, w, lambda, choice_G, Pmat, Bmat) ...
                 + partial_diffusion_SIR(f_tot, x, wk, w, sigma_J, Bmat);  


            % Accumulate the quadrature approximation of \int (C/D) dw

            acc = acc + wj_weights(k) * (Ck ./ Dk); 

        end

        lambda_w(:, i) = acc;
        
    end


end
