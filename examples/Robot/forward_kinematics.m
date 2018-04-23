function [x, y, theta] = forward_kinematics(Ndof, len, qmat, offset)
%FORWARD_KINEMATICS compute forward kinematics
    if nargin < 4
        offset = [0, 0, 0];
    end
    N = size(qmat, 1);
    x = zeros(N, Ndof+1);
    y = zeros(N, Ndof+1);
    theta = zeros(N, Ndof+1);
    % base position
    x(:, 1) = offset(1);
    y(:, 1) = offset(2);
    theta(:, 1) = offset(3);
    for i = 1 : Ndof
        theta(:, i+1) = theta(:, i) + qmat(:, i);
        x(:, i+1) = x(:, i) + len(i) * cos(theta(:, i+1));
        y(:, i+1) = y(:, i) + len(i) * sin(theta(:, i+1));
    end
end

