function d = diffusion_popularity(v_eval, zeta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing diffusion term    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute diffusion term on shifted grid v_eval

d = 0.5 * zeta .* (v_eval.^2);


end