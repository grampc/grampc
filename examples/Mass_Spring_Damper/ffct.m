function dx = ffct(t, x, u, p, userparam)
%#codegen
    m = userparam(1);
    c = userparam(2);
    d = userparam(3);
    
    Nx = size(x, 1);
    NN = Nx/2;
    
    dx = zeros(Nx, 1);
    
    for i=1:NN
        dx(i) = x(NN+i);
    end
    dx(NN+1) = -2*c/m*x(1) + c/m*x(2) - 2*d/m*x(NN+1) + d/m*x(NN+2) + 1/m*u(1);
    
    for i=2:NN-1
        dx(NN+i) = c/m*x(i-1) - 2*c/m*x(i) + c/m*x(i+1) + d/m*x(NN+i-1) - 2*d/m*x(NN+i) + d/m*x(NN+i+1);
    end
    dx(Nx) = c/m*x(NN-1) - 2*c/m*x(NN) + d/m*x(Nx-1) - 2*d/m*x(Nx) - 1/m*u(2);
end