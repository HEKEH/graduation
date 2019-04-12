function f = fb(xb, t, m)
    %omega = 0.4113155;
    f = zeros(length(xb), 1);% * cos(omega * t);
    %f = 10 * 1e-3 * m * [zeros(length(xb) - 10, 1); ones(10, 1)] * cos(omega * t);
end