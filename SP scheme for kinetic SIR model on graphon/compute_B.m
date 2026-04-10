function Bmat = compute_B(x, choice_B)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function computing graphon B   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xcol = x(:);
Nx   = length(xcol);

switch choice_B

    case 1   % B(x,y) = 1

        Bmat = ones(Nx, Nx);

    case 2   % B(x,y) = (x * y)^(-ξ) 

        xi = 0.05;
        cutoff = 1e-5;

        Bmat = (xcol + cutoff).^(-xi) .* (xcol.' + cutoff).^(-xi);

    case 3   % B(x,y) = 1(|x - y| < r) 

        r = 1/5;

        Bmat = abs(xcol - xcol.') < r;

end


end