function px = compute_p(x, choice_p)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Function computing propensity to interact p  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x = x(:);

Nx = length(x);

px = zeros(Nx, 1);

if choice_p == 1

    px(:) = 1;

elseif choice_p == 2   % p(x) = x^(ξ(χ-1)) * ∫_Ω z^(-ξ) * (x^ξ + z^ξ)^(-) dz

    chi = 1/2;

    xi = 0.05;

    cutoff = 1e-10;

    for i = 1:Nx

        xi_val = x(i);

        fx = @(z) (xi_val+cutoff).^(xi*(chi-1)) .* z.^(-xi) .* ((xi_val+cutoff).^xi + z.^xi).^(-chi);

        px(i) = integral(fx, 0, 1, 'ArrayValued', true, 'RelTol', 1e-4, 'AbsTol',1e-8);

    end

elseif choice_p == 3   % p(x) = piecewise expression from appendix A

    chi = 1/2;

    r = 1/5;

    z = linspace(0, 1, 1001).'; dz = z(2) - z(1);

    for i = 1:Nx

        xk = x(i);

        if xk >= 0 && xk < r

            term1 = xk * (1 + (xk + r)/(2*r))^(-chi);
            mask  = (z >= 0) & (z <= r);
            I     = sum( (1 + (xk + r) ./ (z(mask) + r)).^(-chi) ) * dz;
            px(i) = term1 + I;

        elseif xk >= r && xk < 2*r

            term1 = xk * 2^(-chi);
            mask  = (z >= (xk - r)) & (z <= r);
            I     = sum( (1 + (2*r) ./ (z(mask) + r)).^(-chi) ) * dz;
            px(i) = term1 + I;

        elseif xk >= 2*r && xk <= 1 - 2*r

            px(i) = r * 2^(1 - chi);

        elseif xk > 1 - 2*r && xk <= 1 - r

            term1 = (1 - xk) * 2^(-chi);
            mask  = (z >= (1 - r)) & (z <= (xk + r));
            I     = sum( (1 + (2*r) ./ (1 - z(mask) + r)).^(-chi) ) * dz;
            px(i) = term1 + I;

        else

            term1 = (1 - xk) * (1 + (1 - xk + r)/(2*r))^(-chi);
            mask  = (z >= (1 - r)) & (z <= 1);
            I     = sum( (1 + (1 - xk + r) ./ (1 - z(mask) + r)).^(-chi) ) * dz;
            px(i) = term1 + I;

        end

    end

end


end