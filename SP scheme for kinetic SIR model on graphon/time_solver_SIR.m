function f = time_solver_SIR(f, dt, x, w, dw, delta_w, C_w, D_w, scheme_t)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function solving discretized system in time    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[W, ~] = meshgrid(w, x); % Main grid of points (i,j)  --->  size = Nx x Nw

[Nx, Nw] = size(W);


% Compute solution of Fokker-Planck equation

if scheme_t == 1   % Forward Euler

    % Compute solution in intermediate nodes (i+1/2,j)

    f_weighted = (1 - delta_w) .* f(:, 2:Nw) + delta_w .* f(:, 1:Nw-1);


    % Compute numerical fluxes

    Flux_w = C_w .* f_weighted + D_w .* ( f(:, 2:Nw) - f(:, 1:Nw-1) ) / dw;


    % Boundary conditions

    BC = zeros(Nx, 1);


    % Compute solution with FE scheme

    f = f + (dt/dw) * ( [Flux_w, BC] - [BC, Flux_w] );

elseif scheme_t == 2   % Backward Euler

    % Boundary conditions

    BC_w = zeros(Nx, 1);


    % Add boundary conditions to weights delta_w, C_w and D_w

    delta_wBC = [BC_w, delta_w];

    C_wBC = [BC_w, C_w];

    D_wBC = [BC_w, D_w];


    % Transform (Nx x Nw) matrices into (Nx * Nw) vectors

    f_vect = reshape(f', 1, [])';   % size(f) = Nx x Nw  --->  length(f_vect) = Nx * Nw 

    delta_wBC_vect = [reshape(delta_wBC', 1, []), 0]';   % Last zero element needed to complete BC

    C_wBC_vect = [reshape(C_wBC', 1, []), 0]';   % Last zero element needed to complete BC

    D_wBC_vect = [reshape(D_wBC', 1, []), 0]';   % Last zero element needed to complete BC


    % Construct matrix S of linear system

    main_diag = 1 - (dt/dw) * C_wBC_vect(2:end) .* delta_wBC_vect(2:end) ...
                  + (dt/dw) * C_wBC_vect(1:end-1) .* (1 - delta_wBC_vect(1:end-1)) ...
                  + (dt/dw^2) * ( D_wBC_vect(2:end) + D_wBC_vect(1:end-1) );


    upper_diag = - (dt/dw) * C_wBC_vect(2:end-1) .* (1 - delta_wBC_vect(2:end-1)) ...
                 - (dt/dw^2) * D_wBC_vect(2:end-1);

    lower_diag = (dt/dw) * C_wBC_vect(2:end-1) .* delta_wBC_vect(2:end-1) ...
                 - (dt/dw^2) * D_wBC_vect(2:end-1);

    M = diag(main_diag) + diag(upper_diag, +1) + diag(lower_diag, -1);

    S = sparse(M);


    % Solve linear system

    f = S \ f_vect;


    % Transform (Nx * Nw) solution vector into (Nx x Nw) matrix

    f = reshape(f, Nw, Nx)';

end


end