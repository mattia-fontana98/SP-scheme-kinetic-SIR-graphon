function K = drift_SIR(f_tot, x, w_eval, w, lambda, choice_G, Pmat, Bmat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing drift term    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[W_eval, ~] = meshgrid(w_eval, x);   % Shifted grid in w

[Nx, Nw_eval] = size(W_eval);

K = zeros(Nx, Nw_eval);

dy = x(2) - x(1);
dv = w(2) - w(1);

v = w;

for i = 1:Nx

    PBxy = (Pmat(i, :).' .* Bmat(i, :).');   % Pre-computed interaction in x 

    PBf = PBxy .* f_tot;                 

    for j = 1:Nw_eval

        Gwv = compute_G(W_eval(i, j), v, choice_G);


        % Compute drift term on shifted grid w_eval
        
        K(i, j) = lambda * sum(PBf .* (Gwv .* (W_eval(i, j) - v)), 'all') * dv * dy;

    end

end


end