
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Structure-preserving scheme for kinetic SIR model with graphon    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear
clc

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 18);



%% Model parameters


% Discretization in w (opinion of agents)

Nw = 101;                   % Number of nodes
w  = linspace(-1, 1, Nw);   % Domain I = [-1,1]
dw = w(2) - w(1);           % Step size


% Discretization in x (position on graphon)

Nx = 21;                   % Number of nodes
x  = linspace(0, 1, Nx);   % Domain Ω = [0,1]
dx = x(2) - x(1);          % Step size


% Discretization in t (time)

tau = 1;      % Relaxation timescale
Nt  = 2000;   % Number of time steps


% Drift coefficient λ > 0

lambda = 1;


% Diffusion coefficients σ_J^2 > 0

sigma_S = 0.01;
sigma_I = 0.01;
sigma_R = 0.01;


% Epidemiological parameters

beta_contact   = 4/5;   % Baseline contact rate β
gamma_recovery = 3/5;   % Recovery rate γ

R0 = beta_contact / gamma_recovery;   % Basic reproduction number

alpha = 0;   % Opinion-weighted contact exponent



%% Control parameters


% Control plot of total distribution over time

choice_plot = 0;   % 0 --> Deactivate plot
                   %
                   % 1 --> Activate plot


% Type of solver in w (opinion)

scheme_w = 3;   % 1 --> Midpoint (2nd order)
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

choice_IC = 3;   % 1 --> Uniform distribution in x and w
                 %
                 % 2 --> Piecewise constant profile in w
                 %
                 % 3 --> Single rectangular band in w
                 %
                 % 4 --> Uneven blocks, many highly connected positive opinions (Test 3)
                 %
                 % 5 --> Uneven blocks, many poorly connected negative opinions (Test 4)
                

% Type of graphon

choice_B = 1; % 1 --> B(x,y) = 1
              %
              % 2 --> B(x,y) = (x * y)^(-ξ)
              %
              % 3 --> B(x,y) = I(|x - y| < r)
              %
              % Values of ξ ∈ (0,1), and r ∈ (0,1/4) are assigned inside function 'compute_B.m'


% Type of interaction function P

choice_P = 3;   % 1 --> P(x,y) = 1
                %
                % 2 --> P̃(x,y)  = g(p(x)) * g(p(y)) (simplified model, set choice_B = 1)
                %       p(x)    = x^(ξ(χ-1)) * ∫_Ω z^(-ξ) * (x^ξ + z^ξ)^(-) dz
                %       g(p(x)) = p(x) / (a + p(x))
                %
                % 3 --> P̃(x,y)  = g(p(x)) * g(p(y)) (simplified model, set choice_B = 1)
                %       p(x)    = piecewise expression from appendix A
                %       g(p(x)) = p(x) / (a + p(x))
                %
                % 4 --> P(x,y)  = (1 + d_in(x)/d_in(y))^(-χ) (original model, set choice_B = 2)
                %       d_in(x) = ∫_Ω B(x,z) dz
                %       B(x,y)  = (x * y)^(-ξ)
                %
                % 5 --> P(x,y)  = (1 + d(x)/d(y))^(-χ) (original model, set choice_B = 3)
                %       d_in(x) = ∫_Ω B(x,z) dz
                %       B(x,y)  = 1(|x - y| < r)
                %
                % Values of ξ ∈ (0,1), r ∈ (0,1/4), χ > 0, and a > 0 are assigned inside function 'compute_P.m'


% Type of interaction function G

choice_G = 1;   % 1 --> G(w,v) = 1
                %
                % 2 --> G(w,v) = 1(|w - v| < ∆)
                %
                % Value of ∆ ∈ (0,2) is assigned inside function 'define_G.m'


% Option for saving simulation outputs for popularity model

save_data = 0;   % 0 --> Outputs are not saved
                 %
                 % 1 --> Outputs are saved



%% Pre-computing of P and B


% Pre-compute interaction function P

[Pmat, Px_vec] = compute_P(x, choice_P);


% Pre-compute graphon B

Bmat = compute_B(x, choice_B);



%% Define initial data


% Base form of initial distributions

fS_in   = initial_condition_SIR(x, dx, w, dw, choice_IC);
fI_in   = initial_condition_SIR(x, dx, w, dw, 1);
fR_in   = initial_condition_SIR(x, dx, w, dw, 1);


% Macroscopic density of agents in each compartment (rhoS + rhoI + rhoR = 1)

rhoI = 0.001;             % Density of infected agents
rhoR = 0.001;             % Density of recovered agents
rhoS = 1 - rhoI - rhoR;   % Density of susceptible agents

rhoS0 = rhoS;
rhoI0 = rhoI;
rhoR0 = rhoR;

% Initial compartmental distributions

fS = rhoS * fS_in;   % Distribution of susceptible agents
fI = rhoI * fI_in;   % Distribution of infected agents
fR = rhoR * fR_in;   % Distribution of recovered agents

f_tot   = fS + fI + fR;   % Total distribution
ftot_in = f_tot;


% Matrices collecting evolution of compartmental distributions

fS_over_time   = zeros(Nx, Nw, Nt+1);
fI_over_time   = zeros(Nx, Nw, Nt+1);
fR_over_time   = zeros(Nx, Nw, Nt+1);
ftot_over_time = zeros(Nx, Nw, Nt+1);

fS_over_time(:, :, 1)   = fS;   % Register initial distributions
fI_over_time(:, :, 1)   = fI;
fR_over_time(:, :, 1)   = fR;
ftot_over_time(:, :, 1) = ftot_in;


% Vectors collecting evolution of compartmental densitites

rhoS_over_time = zeros(1, Nt+1);
rhoI_over_time = zeros(1, Nt+1);
rhoR_over_time = zeros(1, Nt+1);

rhoS_over_time(1) = rhoS;   % Register initial densities
rhoI_over_time(1) = rhoI;
rhoR_over_time(1) = rhoR;


% Vectors collecting evolution of time and numerical error

time = zeros(1, Nt+1);

if ismember(choice_P, [1, 2, 3]) && choice_G == 1

    % Error committed when approximating steady state

    error = zeros(1, Nt+1);
    error_S = zeros(1, Nt+1);
    error_I = zeros(1, Nt+1);
    error_R = zeros(1, Nt+1);

end


% Vector collecting evolution of effective reproduction number

R_over_time = zeros(1, Nt+1);



%% Compute numerical solution


for j = 1:Nt   % Main loop in time



    % --->  OPINION STEP  <---


    % Introduce distribution and moments of total population

    [rho_tot, mean_tot] = compute_moments_SIR(f_tot, x, dx, w, dw);   % Total density and mean number of agents


    % Compute weights of SP scheme: delta_w, C_w and D_w (in w, pointwise in x)

    [delta_wS, C_wS, D_wS] = compute_weights_SIR(f_tot, x, w, dw, lambda, sigma_S, choice_G, Pmat, Bmat, scheme_w);
    [delta_wI, C_wI, D_wI] = compute_weights_SIR(f_tot, x, w, dw, lambda, sigma_I, choice_G, Pmat, Bmat, scheme_w);
    [delta_wR, C_wR, D_wR] = compute_weights_SIR(f_tot, x, w, dw, lambda, sigma_R, choice_G, Pmat, Bmat, scheme_w);


    % Determine CFL condition for stability and positivity

    if scheme_t == 1   % Stability and positivity of explicit scheme

        dt_S = dw^2 / ( 2 * ( dw * max(max(abs(C_wS))) + max(max(D_wS)) ) );
        dt_I = dw^2 / ( 2 * ( dw * max(max(abs(C_wI))) + max(max(D_wI)) ) );
        dt_R = dw^2 / ( 2 * ( dw * max(max(abs(C_wR))) + max(max(D_wR)) ) );

        dt = min(dt_S, min(dt_I, dt_R));

    elseif scheme_t == 2 % Stability and positivity of IMEX scheme

        dt_S = dw / ( 2 * max(max(abs(C_wS))) );
        dt_I = dw / ( 2 * max(max(abs(C_wI))) );
        dt_R = dw / ( 2 * max(max(abs(C_wR))) );

        dt = min(dt_S, min(dt_I, dt_R));
    end


    % Update solution in each compartment (opinion step)

    fS = time_solver_SIR(fS, dt/tau, x, w, dw, delta_wS, C_wS, D_wS, scheme_t);   % Susceptible agents
    fI = time_solver_SIR(fI, dt/tau, x, w, dw, delta_wI, C_wI, D_wI, scheme_t);   % Infected agents
    fR = time_solver_SIR(fR, dt/tau, x, w, dw, delta_wR, C_wR, D_wR, scheme_t);   % Recovered agents



    % --->  EPIDEMIOLOGICAL STEP  <---


    % Runge-Kutta 4 scheme for SIR

    [RK1_S, RK1_I, RK1_R] = macroscopic_SIR(fS, fI, x, dx, w, dw, beta_contact, gamma_recovery, alpha);
    [RK2_S, RK2_I, RK2_R] = macroscopic_SIR(fS + RK1_S*dt/2, fI + RK1_I*dt/2, x, dx, w, dw, beta_contact, gamma_recovery,alpha);
    [RK3_S, RK3_I, RK3_R] = macroscopic_SIR(fS + RK2_S*dt/2, fI + RK2_I*dt/2, x, dx, w, dw, beta_contact, gamma_recovery,alpha);
    [RK4_S, RK4_I, RK4_R] = macroscopic_SIR(fS + RK3_S*dt, fI + RK3_I*dt, x, dx, w, dw, beta_contact, gamma_recovery,alpha);


    % Update solution in each compartment (epidemiological step)

    fS = fS + dt * (1/6 * RK1_S + 1/3 * RK2_S + 1/3 * RK3_S + 1/6 * RK4_S);   % Susceptible agents
    fI = fI + dt * (1/6 * RK1_I + 1/3 * RK2_I + 1/3 * RK3_I + 1/6 * RK4_I);   % Infected agents
    fR = fR + dt * (1/6 * RK1_R + 1/3 * RK2_R + 1/3 * RK3_R + 1/6 * RK4_R);   % Recovered agents

    f_tot = fS + fI + fR;   % Update total distribution of agents

    
    % Compute moments of compartmental distributions

    [rhoS, mean_wS] = compute_moments_SIR(fS, x, dx, w, dw);
    [rhoI, mean_wI] = compute_moments_SIR(fI, x, dx, w, dw);
    [rhoR, mean_wR] = compute_moments_SIR(fR, x, dx, w, dw);
 
    
    % Record compartmental distributions and moments at this time step

    fS_over_time(:, :, j+1) = fS;
    fI_over_time(:, :, j+1) = fI;
    fR_over_time(:, :, j+1) = fR;
    ftot_over_time(:, :, j+1) = f_tot;

    rhoS_over_time(j+1) = rhoS;
    rhoI_over_time(j+1) = rhoI;
    rhoR_over_time(j+1) = rhoR;


    % Record present time step

    time(j+1) = time(j) + dt;


    % Display conservation of total mass and evolution of total mean

    fprintf('Mass = %.15e  Mean = %f \n step = %f dt =%.3e  \n', ...
        rhoS + rhoI + rhoR, rhoS*mean_wS + rhoI*mean_wI + rhoR*mean_wR, j, dt);


    % Compute steady state of system when known

    if ismember(choice_P, [1, 2, 3]) && choice_G == 1

        [steady_state, rho_P, m_P] = compute_steady_state_SIR(f_tot, ftot_in, x, dx, w, dw, lambda, sigma_S, Px_vec);


        % Compute L1 distance between f_tot and steady state

        error(j+1) = sum(abs(f_tot - steady_state), 'all') * dw * dx;


        % Display conservation of weighted mass and mean in reduced model

        if ismember(choice_P, [2, 3]) && choice_G == 1

            fprintf('rho_P = %.15e  m_P = %f \n ', rho_P, m_P);

        end

    end


    % Plot total distribution over time

    if choice_plot == 1

        figure(1)
        plot(w(2:end-1), f_tot(1, 2:end-1),'o-b');
        hold on

        if ismember(choice_P, [1, 2, 3]) && choice_G == 1

            plot(w, steady_state(1, :), 'o-r');

        end

        hold off
        xlabel('v');
        ylabel('f(v)');
        title(sprintf('Density = %f Time = %f', rhoS+rhoI+rhoR, (j+1) * dt));

        drawnow

    end


    % Compute effective reproduction number

    if alpha == 1

        Salpha = sum(fS * (1 - w(:)).^alpha, 'all') * dx * dw;    % = ∫_(Ω x I) (1-w)^α fS(t,x,w) dx dw = ρS(t) * (1-mS(t))
        Malpha = sum(fI * (1 - w(:)).^alpha, 'all') * dx * dw;    % = ∫_(Ω x I) (1-w)^α fI(t,x,w) dx dw = ρI(t) * (1-mI(t))

        R_eff = beta_contact * Salpha * Malpha / (gamma_recovery * rhoI);   % = β / γ * (1-mS(t)) * (1-mI(t)) * ρS(t)

        R_over_time(j+1) = R_eff;

    end


end   % End of main loop in time



%% Convergence to steady state

if ismember(choice_P, [1, 2, 3]) && choice_G == 1

    % Compute final density of susceptibles

    f = @(S) log(S/rhoS0) + R0*(rhoS0 + rhoI0 - S);

    b = rhoS0 * (1 - 1e-12);

    rhoS_inf = fzero(f, [1e-12, b]);
    rhoI_inf = 0;
    rhoR_inf = 1 - rhoS_inf;


    % Compute steady states

    [B_S, ~, ~] = compute_steady_state_SIR(f_tot, fS_in, x, dx, w, dw, lambda, sigma_S, Px_vec);
    [B_R, ~, ~] = compute_steady_state_SIR(f_tot, fR_in, x, dx, w, dw, lambda, sigma_R, Px_vec);

    steady_state_S = rhoS_inf * B_S;
    steady_state_I = zeros(size(fI));
    steady_state_R = rhoR_inf * B_R;


    % Compute errors

    iw = 2:length(w)-1;   % Use 1:length(w) to include boundaries

    for k = 1:Nt+1

        error_S(k) = sum(abs(fS_over_time(:, iw, k) - steady_state_S(:, iw)), 'all') * dw * dx;
        error_I(k) = sum(abs(fI_over_time(:, iw, k) - steady_state_I(:, iw)), 'all') * dw * dx;
        error_R(k) = sum(abs(fR_over_time(:, iw, k) - steady_state_R(:, iw)), 'all') * dw * dx;

    end

end



%% Save outputs in .mat files


% Save data to run simulations of popularity model, mandatory for Test 2

if save_data == 1

    if choice_IC == 2   % Test 1: Piecewise constant profile in w

        save Output_kineticSIR_piecewise_w.mat x w ftot_over_time time

    elseif choice_IC == 3   % Test 1: Single rectangular band in w

        save Output_kineticSIR_band_w.mat x w ftot_over_time time

    end

end


% Automatically save whole workspace, mandatory for plotting

save Output_kineticSIR.mat
