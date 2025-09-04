function dx = ffct(t, x, u, p, userparam)
%#codegen
    pSys = userparam;
    
    dx = zeros(6, 1);
    dx(1) = x(2);
    dx(2) = -(sin(x(5))*u(1)) + cos(x(5))*pSys(1)*u(2);
    dx(3) = x(4);
    dx(4) = -1 + cos(x(5))*u(1) + pSys(1)*sin(x(5))*u(2);
    dx(5) = x(6);
    dx(6) = u(2);
end