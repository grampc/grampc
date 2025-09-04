function dx = ffct(t, x, u, p, userparam)
%#codegen
    dx = zeros(6, 1);
    dx(1) = x(2);
    dx(2) = u(1);
    dx(3) = x(4);
    dx(4) = u(2);
    dx(5) = x(6);
    dx(6) = -((9.81 * sin(x(5)) + cos(x(5))*u(1) + 2 * x(4) * x(6)) / x(3));
end