function dx = ffct(t, x, u, p, userparam)
%#codegen
    pSys = userparam;
    
    dx = zeros(5, 1);
    dx(1) = cos(x(3))*x(5);
    dx(2) = sin(x(3))*x(5);
    dx(3) = (x(4)*x(5))/(pSys(1)*(1 + 1/pSys(2)^2*x(5)^2));
    dx(4) = u(1);
    dx(5) = u(2);
end