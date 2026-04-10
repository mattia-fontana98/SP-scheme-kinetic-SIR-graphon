function G = compute_G(w_mid, v, choice_G)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing interaction function G    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nv = length(v);

if choice_G == 1   % G(w,v)= 1

    G = ones(1, Nv);
    
elseif choice_G == 2   % G(w,v) = 1(|w - v| < ∆) 

    Delta = 0.5;

    G = abs(w_mid - v) < Delta;

end


end