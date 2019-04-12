function f = fc(xc, t)
    omega = 2;
    f = 1e-3 * sin(pi * xc)' * cos(omega * t);
    %f = zeros(length(xc), 1);
end