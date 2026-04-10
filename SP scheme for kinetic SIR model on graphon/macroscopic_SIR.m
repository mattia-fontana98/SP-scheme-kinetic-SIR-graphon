function [RK_S, RK_I, RK_R] = macroscopic_SIR(f_S, f_I, x, dx, w, dw, beta_contact, gamma_recovery, alpha)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing Runge-Kutta weights for SIR model    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[W, ~] = meshgrid(w, x);

M_alpha = sum(((1 - W).^alpha) .* f_I, 'all') * dw * dx;

contact_term = beta_contact * ((1 - W).^alpha) .* f_S .* M_alpha;

recovery_term = gamma_recovery * f_I;

RK_S = -contact_term;
RK_I = contact_term - recovery_term;
RK_R = recovery_term;


end