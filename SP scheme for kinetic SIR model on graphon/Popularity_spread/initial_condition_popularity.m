function h0 = initial_condition_popularity(v, dv, L, choice_IC)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function defining initial distribution     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h0 = zeros(size(v));
v0 = L/3;

switch choice_IC

    case 1   % Left step distribution

        h0(v <= v0) = 1;

    case 2   % Decreasing distribution on [0,L/3]

        mask = (v >= 0) & (v <= v0);
        h0(mask) = 1 - v(mask)/v0;

end


% Normalization

rho0 = sum(h0) * dv;

h0 = h0 / rho0;


end