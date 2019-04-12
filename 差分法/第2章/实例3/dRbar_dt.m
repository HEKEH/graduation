function ret = dRbar_dt(t, R_bar, R_bar_p, Nc, Nb, K_of_fR, row_vector_fR6, H, lambda, Lc, Bc, dxc, theta, Kc1, Kc2, Kc3, Kb1, Kb2, Fc, Fb)
    [Rc, Rb] = fR(R_bar, Nc, Nb, K_of_fR, row_vector_fR6);
    wcNc = Rc(end);
    WcNc = ones(Nc - 1, 1) * wcNc;
    GAMMA_c = get_GAMMA_c(Rc, H, lambda, Lc, Bc, dxc, theta, WcNc, Nc);
    dRcbar_dt = Kc1 * Rc + Kc2 * Fc(t) + Kc3 * WcNc + Kc2 * GAMMA_c;
    dRbbar_dt = Kb1 * Rb + Kb2 * Fb(t);
    ret = R_bar_p - [dRcbar_dt; dRbbar_dt];
end