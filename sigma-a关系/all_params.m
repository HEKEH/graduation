Configs;
Ib = K * Ec * Ac^2 / Eb; %m4
beta1 = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);
k1 =  m / (beta1 * cos(theta));
k2 =  m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta));
k3 =  m * (- 2 * H * tan(theta) + 1)^2 * lambda^2 * cos(theta) / (4 * beta1);
k4 =  -m / (beta1 * cos(theta));
phi_params;
psi2_params;
psi3_params;
int_func = @(idx, func)integral(@(xc)func(xc, idx), 0, 1);
PHI_c31 = @(xc, idx)H * lambda^2 * (phi_c_p2(xc, idx) * int_func(idx, psi_c2) + psi_c2_p2(xc, idx) * int_func(idx, phi_c));
PHI_c32 = @(xc, idx)H * lambda^2 * (phi_c_p2(xc, idx) * int_func(idx, psi_c3) + psi_c3_p2(xc, idx) * int_func(idx, phi_c));
PHI_c33 = @(xc, idx)H * lambda^2 * (1/2 - H * tan(theta)) * (phi_c_p2(xc, idx) .* psi_c2(1, idx) + psi_c2_p2(xc, idx) .* phi_c(1, idx));
PHI_c34 = @(xc, idx)H * lambda^2 * (1/2 - H * tan(theta)) * (phi_c_p2(xc, idx) .* psi_c3(1, idx) + psi_c3_p2(xc, idx) .* phi_c(1, idx));
I_c31 = @(idx)2 * H * lambda^2 *int_func(idx, @ (xc, idx)phi_c_p(xc, idx) .* psi_c2_p (xc, idx));
I_c32 = @(idx)2 * H * lambda^2 *int_func(idx, @ (xc, idx)phi_c_p(xc, idx) .* psi_c3_p (xc, idx));
PHI_c3 = @(xc, idx) PHI_c31(xc, idx) - PHI_c33(xc, idx) - I_c31(idx) + PHI_c32(xc, idx) - PHI_c34(xc, idx) - I_c32(idx);

clear i