function [Nv, L] = select_Nv_L(alpha, beta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function defining number of grid points Nv and domain length L    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


npp         = 20;      % Number of points used to resolve peak
epsTail     = 1e-12;   % Tolerance on tail mass
Nv_min      = 101;     % Minimum number of grid points
Nv_max      = 801;     % Maximum number of grid points
Lmax        = 1e6;     % mMximum allowed domain length
Lmin_factor = 8;       % Minimum domain size as a multiple of peak location


% Approximate width of peak

v_peak     = beta/(alpha + 1);
sigma_peak = beta/(alpha + 1)^(3/2);


% Target mesh size for suitable resolution of peak

dv_target = sigma_peak / npp;


% Minimum domain size

L_min = Lmin_factor * v_peak;


% Domain size required to control tail mass

y_eps  = gammaincinv(epsTail, alpha, 'lower');
L_tail = beta / y_eps;


% Initial domain length after imposing minimum and maximum bounds

L0 = min(max(L_tail, L_min), Lmax);


% Maximum admissible domain length compatible with dv <= dv_target, under constraint Nv <= Nv_max

L_target = min(Lmax, (Nv_max - 1) * dv_target);

if L_min > L_target

    Nv = Nv_max;
    L  = min(L_min, Lmax);
    
    return;

end


% Final domain length

L = min(L0, L_target);


% Required number of grid points to achieve target mesh size

Nv_need = ceil(L / dv_target) + 1;
Nv = min(max(Nv_need, Nv_min), Nv_max);


end