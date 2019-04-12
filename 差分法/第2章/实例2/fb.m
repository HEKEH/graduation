function f = fb(xb, t, m)
    omega = 2;
    f = 1e-4 * sin(pi * xb)' * cos(omega * t);
    %f = 10 * 1e-3 * m * [zeros(length(xb) - 10, 1); ones(10, 1)] * cos(omega * t);
end