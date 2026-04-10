
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Structure-preserving scheme for popularity model    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear
clc

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 18);


% Load necessary output data from kinetic SIR simulations

load(fullfile('..', 'Output_kineticSIR_piecewise_w.mat'));



%% Model parameters


% Discretization in w (opinion of agents)

Nw = length(w);
dw = w(2) - w(1);   % Step size


% Discretization in x (position on graphon)

Nx = length(x);
dx = x(2) - x(1);   % Step size


% Discretization in t (time)

tau = 1; % Relaxation timescale


% Drift coefficient μ > 0

mu = 1.5;


% Diffusion coefficient ζ^2 > 0

zeta = 1;


% Spread factor θ > 0

theta = 5;


% Hypothesized danger level of disease w_hat ∈ [−1,1]

w_hat = 0.3;



%% Control parameters


% Type of solver in v (popularity)

scheme_v = 1;   % 1 --> Midpoint (2nd order)
                %
                % 2 --> Milne (4th order)
                %
                % 3 --> Newton-Cotes (6th order)
                %
                % 4 --> Gauss-Legendre (8th order)


% Type of solver in t (time)

scheme_t = 2;   % 1 --> Forward Euler (explicit, not AP w.r.t. tau)
                %
                % 2 --> Backward Euler (IMEX, AP w.r.t. tau)


% Type of initial condition

choice_IC = 1;   % 1 --> Left step distribution
                 %
                 % 2 --> Decreasing distribution on [0,L/3]


% Type of propensity to interact p

choice_p = 3;   % 1 --> p(x) = 1
                %
                % 2 --> p(x) = x^(ξ(χ-1)) * int_Ω z^(-ξ) * (x^ξ + z^ξ)^(-χ) dz
                %
                % 3 --> p(x) = piecewise expression from Appendix A
                %
                % Values of ξ ∈ (0,1), r ∈ (0,1/4), and χ > 0 are assigned inside function 'compute_p.m'



%% Pre-computing of p and F


% Pre-compute propensity to interact p

p_x = compute_p(x, choice_p);

Nt_SIR = size(ftot_over_time, 3);


% Pre-compute nonlocal term F

F_series = zeros(1, Nt_SIR);

for k = 1:Nt_SIR

    F_series(k) = compute_F(ftot_over_time(:, :, k), w, dx, dw, p_x, w_hat);

end


% Compute F_inf

F_inf = compute_F(ftot_over_time(:, :, end), w, dx, dw, p_x, w_hat);



%% Discretization in v (popularity of a product)


% Compute inverse gamma parameters

alpha = 1 + 2*mu/zeta;
beta  = 2 * theta * F_inf/zeta;


% Define number of nodes Nv and domain upper length L

[Nv, L] = select_Nv_L(alpha, beta);

v  = linspace(0, L, Nv);   % Domain [0,L]
dv = v(2) - v(1);          % Step size



%% Define initial data


h0 = initial_condition_popularity(v, dv, L, choice_IC);

h_save = zeros(Nv, Nt_SIR);
h_save(:, 1) = h0;

tH = 0;
h  = h0;

t_final = time(end);

kSave = 2;



%% Compute numerical solution


while tH < t_final - 1e-14   % Main loop in time


    F_now = F_series(max(1, kSave-1));


    % Compute weights of SP scheme: delta_v, C_v and D_v

    [delta_v, C_v, D_v] = compute_weights_popularity(F_now, v, dv, mu, zeta, theta, scheme_v);


    % Determine CFL condition for stability and positivity

    C = max(max(abs(C_v)));

    D = max(max(D_v));


    if scheme_t == 1 % Stability and positivity of explicit scheme

        dt = dv^2 / (2 * (C * dv + D));

    elseif scheme_t == 2 % Stability and positivity of IMEX scheme

        dt = dv / (2 * C);

    end


    tolT = 1e-14;


    % Synchronize popularity time-stepping with kinetic SIR output times

    if kSave <= Nt_SIR

        dt_event = time(kSave) - tH;

        if dt_event <= tolT

            h_save(:, kSave) = h;
            kSave = kSave + 1;

            continue
        end

        dt = min(dt, dt_event);
    end


    dt = min(dt, t_final - tH);


    % Compute solution

    h  = time_solver_popularity(h, dt/tau, v, dv, delta_v, C_v, D_v, scheme_t);


    % Update time

    tH = tH + dt;


    % Check next SIR snapshot after this step is reached

    if kSave <= Nt_SIR && (tH >= time(kSave) - tolT)

        h_save(:, kSave) = h;
        kSave = kSave + 1;

    end


    % Compute moments of solution

    [rho, mean_v] = compute_moments_popularity(h, v, dv);


    % Display conservation of total mass and evolution of total mean

    fprintf('Total mass = %.15e  Mean = %f \n ', rho, mean_v);


end   % End of main loop in time



%% Clear command window


clc



%% Verify convergence to steady state


h_steady = compute_steady_state_popularity(v, dv, mu, zeta, theta, F_inf);

tH = time(:);

h_over_time = h_save;

F_used = F_series(:);

errL1 = sum(abs(h_over_time -  h_steady(:)), 1)' * dv;

eF = abs(F_series(:) - F_inf);



%% Save data for plot

save Output_popularity.mat
