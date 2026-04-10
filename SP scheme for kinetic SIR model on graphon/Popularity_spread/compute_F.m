function F = compute_F(f_tot, w, dx, dw, p_x, w_hat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Function computing nonlocal term F   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute F(t) = ∫_(Ω x I) p(x) f_tot(t,x,w) 1(w >= w_hat) dx dw 

Nw = length(w);

iw = 2:Nw-1; % Use 1:Nw to include boundaries 

H = (w(iw) >= w_hat);

F = sum((p_x(:) .* ones(1, length(iw))) .* (f_tot(:, iw) .* H), 'all') * dw * dx;


end