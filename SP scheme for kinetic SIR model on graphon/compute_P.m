function [Pmat, Px_vec] = compute_P(x, choice_P)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Function computing interaction function P     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xcol = x(:);
Nx   = length(xcol);

Pmat   = [];
Px_vec = [];

switch choice_P

    case 1   % P(x,y) = 1

        Pmat   = ones(Nx, Nx);
        Px_vec = ones(Nx, 1);
       

    case 2   % P̃(x,y)  = g(p(x)) * g(p(y))
             % p(x)    = x^(ξ(χ-1)) * ∫_Ω {z^(-ξ) * (x^ξ + z^ξ)^(-χ)} dz
             % g(p(x)) = p(x) / (a + p(x))

        xi  = 0.25;
        chi = 2;
        a   = 1;
        cutoff = 1e-10;

        fx = @(z) (xcol + cutoff).^(xi*(chi - 1)) .* ...
            z.^(-xi) .* ((xcol + cutoff).^xi + z.^xi).^(-chi);

        px  = integral(fx, 0, 1, 'ArrayValued', true, 'RelTol', 1e-4, 'AbsTol', 1e-8);

        gpx = px ./ (a + px);

        Pmat = gpx * (gpx.');

        Px_vec = gpx;


    case 3   % P̃(x,y)  = g(p(x)) * g(p(y))
             % p(x)    = piecewise expression from Appendix A
             % g(p(x)) = p(x) / (a + p(x))

        r   = 1/5;
        chi = 1/2;
        a   = 1/2;

        z  = xcol;
        dz = x(2) - x(1);

        px = zeros(Nx, 1);

        for ix = 1:Nx
            xk = xcol(ix);

            if xk >= 0 && xk < r

                term1 = xk * (1 + (xk + r)/(2*r))^(-chi);
                mask  = (z >= 0) & (z <= r);
                I     = sum( (1 + (xk + r) ./ (z(mask) + r)).^(-chi) ) * dz;
                px(ix) = term1 + I;

            elseif xk >= r && xk < 2*r

                term1 = xk * 2^(-chi);
                mask  = (z >= (xk - r)) & (z <= r);
                I     = sum( (1 + (2*r) ./ (z(mask) + r)).^(-chi) ) * dz;
                px(ix) = term1 + I;

            elseif xk >= 2*r && xk <= 1 - 2*r

                px(ix) = r * 2^(1 - chi);

            elseif xk > 1 - 2*r && xk <= 1 - r

                term1 = (1 - xk) * 2^(-chi);
                mask  = (z >= (1 - r)) & (z <= (xk + r));
                I     = sum( (1 + (2*r) ./ (1 - z(mask) + r)).^(-chi) ) * dz;
                px(ix) = term1 + I;

            else

                term1 = (1 - xk) * (1 + (1 - xk + r)/(2*r))^(-chi);
                mask  = (z >= (1 - r)) & (z <= 1);
                I     = sum( (1 + (1 - xk + r) ./ (1 - z(mask) + r)).^(-chi) ) * dz;
                px(ix) = term1 + I;

            end

        end

        gpx = px ./ (a + px);

        Pmat = gpx * (gpx.');

        Px_vec = gpx;


    case 4   % P(x,y)  = (1 + d_in(x)/d_in(y))^(-χ)
             % d_in(x) = ∫_Ω B(x,z) dz
             % B(x,z)  = (x * z)^(-ξ)

        xi  = 0.05;
        chi = 2;
        cutoff = 1e-5;

        d_in = zeros(Nx, 1);

        for ix = 1:Nx

            xk = xcol(ix);

            Bxz = @(z) (xk + cutoff*(xk==0)).^(-xi) .* z.^(-xi);

            d_in(ix) = integral(Bxz, 0, 1);

        end

        Pmat = (1 + d_in ./ (d_in.')).^(-chi);
        Px_vec = zeros(Nx, 1);


    case 5   % P(x,y)  = (1 + d_in(x)/d_in(y))^(-χ)
             % d_in(x) = ∫_Ω B(x,z) dz
             % B(x,z)  = 1(|x - z| < r)

        r   = 1/5;
        chi = 2;

        d_in = zeros(Nx, 1);

        for ix = 1:Nx

            xk = xcol(ix);
            Bxz = @(z) double(abs(xk - z) < r);
            d_in(ix) = integral(Bxz, 0, 1);
            
        end

        Pmat   = (1 + d_in ./ (d_in.')).^(-chi);
        Px_vec = zeros(Nx, 1);
        
end


end