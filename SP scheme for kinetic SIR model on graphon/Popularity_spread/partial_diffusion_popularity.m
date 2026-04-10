function d = partial_diffusion_popularity(v_eval, zeta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing partial derivative of diffusion term    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute partial derivative of diffusion term on shifted grid v_eval

d = zeta .* v_eval;


end