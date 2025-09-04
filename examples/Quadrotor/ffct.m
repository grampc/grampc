function dx = ffct(t, x, u, p, userparam)
%#codegen
    dx = zeros(9, 1);
    dx(1) = x(2);
    dx(2) = u(1) * (cos(x(7))*cos(x(9))*sin(x(8)) + sin(x(9))*sin(x(7)));
    dx(3) = x(4);
    dx(4) = u(1) * (sin(x(9))*cos(x(7))*sin(x(8)) - cos(x(9))*sin(x(7)));
    dx(5) = x(6);
    dx(6) = -9.81 + u(1)*cos(x(8))*cos(x(7));
    dx(7) = (u(2)*cos(x(7)) + u(3)*sin(x(7)))/cos(x(8));
    dx(8) = u(3)*cos(x(7)) - u(2)*sin(x(7));
    dx(9) = u(2)*cos(x(7))*tan(x(8)) + u(3)*tan(x(8))*sin(x(7)) + u(4);
end