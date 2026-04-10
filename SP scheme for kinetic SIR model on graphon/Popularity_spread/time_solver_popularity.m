function f = time_solver_popularity(f, dt, v, dv, delta_v, C_v, D_v, scheme_t)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function solving discretized system in time    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nv = length(v);


% Compute solution of Fokker-Planck equation

if scheme_t == 1   % Forward Euler

    % Compute solution in intermediate nodes k+1/2

    f_weighted =  (1 - delta_v) .* f(2:Nv) + delta_v .* f(1:Nv-1);


    % Compute numerical fluxes

    Flux_v = C_v .* f_weighted + D_v .* (f(2:Nv) - f(1:Nv-1)) / dv;


    % Boundary conditions

    BC = 0;


    % Compute solution with FE scheme

    f = f + (dt/dv) * ( [Flux_v, BC] - [BC, Flux_v] );

elseif scheme_t == 2   % Backward Euler

    % Boundary conditions

    BC = 0;


    % Add boundary conditions to weights delta_v, C_v and D_v

    delta_vBC = [BC, delta_v];

    C_vBC = [BC, C_v];

    D_vBC = [BC, D_v];

    f_vect = f';

    delta_vBC_vect = [delta_vBC, 0]';   % Last zero element needed to complete BC

    C_vBC_vect = [C_vBC, 0]';   % Last zero element needed to complete BC

    D_vBC_vect = [D_vBC, 0]';   % Last zero element needed to complete BC


    % Construct matrix S of linear system

    main_diag = 1 - (dt/dv) * C_vBC_vect(2:end) .* delta_vBC_vect(2:end) ...
                  + (dt/dv) * C_vBC_vect(1:end-1) .* (1 - delta_vBC_vect(1:end-1)) ...
                  + (dt/dv^2) * ( D_vBC_vect(2:end) + D_vBC_vect(1:end-1) );


    upper_diag = - (dt/dv) * C_vBC_vect(2:end-1) .* (1 - delta_vBC_vect(2:end-1)) ...
                 - (dt/dv^2) * D_vBC_vect(2:end-1);

    lower_diag = (dt/dv) * C_vBC_vect(2:end-1) .* delta_vBC_vect(2:end-1) ...
                 - (dt/dv^2) * D_vBC_vect(2:end-1);


    M = diag(main_diag) + diag(upper_diag, +1) + diag(lower_diag, -1);

    S = sparse(M);


    % Solve linear system

    f = S \ f_vect;


    % Transform solution column vector into row vector

    f = f';

end


end