clear;
func_sym;
syms lambda H

params;
psi_ck_sym = psi_ck_known_sym + double(eval(C_k_sym))' * psi_ck_uncertain_sym';
%temp_ret = -omega_m^2 * psi_ck_known_sym - diff(psi_ck_known_sym, xc, 2) + lambda^2 * int(psi_ck_known_sym, xc, 0, 1) + lambda^2 * (H*tan(theta) - 1/2) * subs(psi_ck_known_sym, xc, 1);
temp_ret = -omega_m^2 * psi_ck_sym - diff(psi_ck_sym, xc, 2) + lambda^2 * int(psi_ck_sym, xc, 0, 1) + lambda^2 * (H*tan(theta) - 1/2) * subs(psi_ck_sym, xc, 1);

psi_ck = matlabFunction(eval(psi_ck_sym));
temp_ret_func = matlabFunction(eval(temp_ret));
temp_func =@(xc) temp_ret_func(xc) - PHI2(xc)+PHI5(xc)+I2 + m * J_m * phi_c(xc, m_idx);
temp_func(0.5)
psi_ck(1)