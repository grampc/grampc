function plot_configuration(Ndof, len, q, offset, c)
%PLOT_CONFIGURATION plot configuration (hold on!)
    if nargin < 5
        c = 'b';
    end
    
    [x, y, ~] = forward_kinematics(Ndof, len, q, offset);
    plot(x(1), y(1), 'o', 'Color', c);
    for i = 1 : length(q)
        line([x(i) x(i+1)], [y(i) y(i+1)], 'LineWidth', 2, 'Color', c);
        plot(x(i+1), y(i+1), 'o', 'Color', c);
    end
end

