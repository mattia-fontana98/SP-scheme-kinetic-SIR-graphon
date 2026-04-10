function z = drift_popularity(F, v_eval, mu, theta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      Function computing drift term      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute drift term on shifted grid v_eval

z = mu .* v_eval - theta .* F;

end
