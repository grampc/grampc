function dx = ffct(t, x, u, p, userparam)
%#codegen
    pSys = userparam;
    dx = zeros(4, 1);
    
    k1 = pSys(1) * exp(-(pSys(3)/(273.15 + x(3))));
    k2 = pSys(2) * exp(-(pSys(4)/(273.15 + x(3))));
    h = -pSys(8) * (k1*pSys(11)*x(1) + k1*pSys(12)*x(2) + k2*x(1)^2*pSys(13));
    
    dx(1) = (-(k2*x(1)^2) + u(1)*(pSys(10) - x(1)) - k1*x(1))/pSys(14);
    dx(2) = (k1*x(1) - k1*x(2) - u(1)*x(2))/pSys(14);
    dx(3) = (h + u(1)*(pSys(9) - x(3)) + pSys(5)*(-x(3) + x(4)))/pSys(14);
    dx(4) = (pSys(7)*u(2) + pSys(6)*(x(3)-x(4)))/pSys(14);
end