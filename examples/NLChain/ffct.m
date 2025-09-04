function dx = ffct(t, x, u, p, userparam)
%#codegen
    N = 4;
    g = [0.0; 0.0; 9.81];
    m = 0.45;
    Nx = 3*(2*N-1);
    dx = zeros(Nx, 1);

    % derivatives of x
    dx(1:3*(N-1)) = x(3*N+1:end);
    dx(3*(N-1)+1:3*N) = u;

    % derivatives of v
    dv = zeros((N-1)*3, 1);
    for i=1:(N-1)
        if i == 1
            dv((i-1)*3+1:(i-1)*3+3) = (F(x((i)*3+1:(i)*3+3), x((i-1)*3+1:(i-1)*3+3), N) - F(x((i-1)*3+1:(i-1)*3+3), [0.0; 0.0; 0.0], N)) * N/m - g;
        else
            dv((i-1)*3+1:(i-1)*3+3) = (F(x((i)*3+1:(i)*3+3), x((i-1)*3+1:(i-1)*3+3), N) - F(x((i-1)*3+1:(i-1)*3+3), x((i-2)*3+1:(i-2)*3+3), N)) * N/m - g;
        end
    end
    dx(3*N+1:end) = dv;
end

function out = F(x1, x0, N)
    k = 0.1;
    lr = 0.55;
    out = (x1-x0) * k * (N - lr/norm(x1-x0));
end