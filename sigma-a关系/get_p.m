function p = get_p(omegas_num, omegas_idx, m, theta, PHI_c3, phi_c, phi_b, Fc, Fb)
    p = zeros(omegas_num, 4);

    for idx = omegas_idx
        p(idx, 1) = m * integral(@(xc)phi_c(xc, idx).^ 2, 0, 1);
        p(idx, 2) = m * integral(@(xc)phi_c(xc, idx) .* PHI_c3(xc, idx), 0, 1);
        p(idx, 3) = 1/2 * m * integral(@(xc)phi_c(xc, idx).* Fc(xc), 0, 1);
        p(idx, 4) = 1/2 * cos(theta)^3 * integral(@(xb)phi_b(xb, idx).* Fb(xb), 0, 1);
    end
end


