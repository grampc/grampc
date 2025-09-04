function dx = ffct(t, x, u, p, userparam)
%#codegen
    pSys = userparam;
    
    dx = zeros(4, 1);
    dx(1) = (u(1) - pSys(1)*x(1) + pSys(3)*x(2)*x(3))/pSys(2);
    dx(2) = (u(2) - pSys(1)*x(2) - (pSys(4) + pSys(2)*x(1))*x(3))/pSys(3);
    dx(3) = (0.5*(-2*pSys(5) * pSys(8) + 3 * pSys(5)^2*(pSys(4) + (pSys(2) - pSys(3))*x(1))*x(2) - 2 * pSys(7) * x(3) )) / pSys(6);
    dx(4) = x(3);
end