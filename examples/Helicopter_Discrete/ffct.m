function dx = ffct(t, x, u, p, userparam)
%#codegen
    dx = zeros(6, 1);
    dx(1) = x(3);
    dx(2) = 0.01*u(1) - 0.01*u(2) + x(4);
    dx(3) = 0.19*u(1) + 0.19*u(2);
    dx(4) = 1.32*u(1) - 1.32*u(2);
    dx(5) = x(1);
    dx(6) = x(2);
end