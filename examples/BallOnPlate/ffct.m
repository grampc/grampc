function dx = ffct(t, x, u, p, userparam)
%#codegen
    dx = zeros(2, 1);
    dx(1) = -0.04*u(1) + x(2);
    dx(2) = -7.01*u(1);
end