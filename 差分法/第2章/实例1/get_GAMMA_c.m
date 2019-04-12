function GAMMA_c = get_GAMMA_c(Rc, H, lambda, Lc, Bc, dxc, theta, WcNc, Nc)
    Wc = Rc(1: Nc + 1, :);
    temp = (Wc(3: end) - Wc(1: end - 2));
    sqr_sum = temp' * temp;
    Xc = ones(Nc - 1, 1) * sqr_sum;
    k = H * lambda^2;
    GAMMA_c = k * (Lc * Wc) .* (Bc * Wc) - k / 4 / dxc * Xc + k * (H * tan(theta) - 1/2) * (Lc * Wc) .* WcNc;
end