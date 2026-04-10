function f_in = initial_condition_SIR(x, dx, w, dw, ic)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function defining initial distributions    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[W, X] = meshgrid(w, x);
[Nx, Nw] = size(X);

if ic == 1   % Uniform distribution in x and w

    f_in = 1/2 * ones(Nx, Nw);

elseif ic == 2   % Piecewise constant profile in w

    levels = [2, 1];

    g = zeros(1, Nw);
    g(w < 0)  = levels(1);
    g(w >= 0) = levels(2);

    f_in = repmat(g, Nx, 1);

elseif ic == 3   % Single rectangular band in w

    f_in = zeros(Nx, Nw);

    l = 1/3;
    c = 1.0;
    center = 0.5;

    band = (W >= center - l/2) & (W <= center + l/2);
    f_in(band) = c;


elseif ic == 4   % Uneven blocks, many highly connected positive opinions (Test 3)

    f_in = zeros(Nx, Nw);

    ax1 = 0.2;   aw1 = 0.25;
    x0_1 = 0.8;  w0_1 = -0.6;
    mask1 = (abs(X - x0_1) <= ax1) & (abs(W - w0_1) <= aw1);

    ax2 = 0.2;   aw2 = 0.25;
    x0_2 = 0.2;  w0_2 = 0.6;
    mask2 = (abs(X - x0_2) <= ax2) & (abs(W - w0_2) <= aw2);

    f_in(mask1) = 1/4;
    f_in(mask2) = 1.0;

elseif ic == 5   % Uneven blocks, many poorly connected negative opinions (Test 4)

    f_in = zeros(Nx, Nw);

    ax1 = 0.2;   aw1 = 0.3;
    x0_1 = 0.8;  w0_1 = -0.6;
    mask1 = (abs(X - x0_1) <= ax1) & (abs(W - w0_1) <= aw1);

    ax2 = 0.15;   aw2 = 0.2;
    x0_2 = 0.15;  w0_2 = 0.7;
    mask2 = (abs(X - x0_2) <= ax2) & (abs(W - w0_2) <= aw2);

    f_in(mask1) = 1.0;
    f_in(mask2) = 1/4;

end


% Normalization

rho_in = sum(f_in, 'all') * dw * dx;

f_in = f_in / rho_in;


end