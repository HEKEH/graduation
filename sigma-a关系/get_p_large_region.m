p = get_p(omegas_num, omegas_idx, m, theta, PHI_c3, phi_c, phi_b, Fc, Fb);
OMEGA_range = [];
start = 0;
for idx = 1: omegas_num
    mid = ceil(omegas(idx) / 1e-3) * 1e-3;
    section = [start: 1e-3: mid - 5e-3, mid - 5e-3 + 1e-5: 1e-5: mid + 5e-3 - 1e-5];
    OMEGA_range = [OMEGA_range, section];
    start = mid + 5e-3;
end
section = mid + 5e-3 : 1e-3: mid + 200e-3;
OMEGA_range = [OMEGA_range, section];
sigma_range = zeros(omegas_num, length(OMEGA_range));
for idx = 1: omegas_num
    sigma_range(idx, :) = OMEGA_range - omegas(idx);
end
sigma_len = size(sigma_range, 2);



