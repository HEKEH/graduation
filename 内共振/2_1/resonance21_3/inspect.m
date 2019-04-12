% clear;
% func_sym;
% syms lambda H
% params;
% C_k = double(eval(C_k_sym))';
% %C_k = [1e13,-5e15,1,1,1,1];
% psi_ck_sym = psi_ck_known_sym + C_k * psi_ck_uncertain_sym';
% %temp_ret = -omega_n^2 * psi_ck_known_sym - diff(psi_ck_known_sym, xc, 2) + lambda^2 * int(psi_ck_known_sym, xc, 0, 1) + lambda^2 * (H*tan(theta) - 1/2) * subs(psi_ck_known_sym, xc, 1);
% temp_ret = -omega_n^2 * psi_ck_sym - diff(psi_ck_sym, xc, 2) + lambda^2 * int(psi_ck_sym, xc, 0, 1) + lambda^2 * (H*tan(theta) - 1/2) * subs(psi_ck_sym, xc, 1);
% 
% temp_ret_func = matlabFunction(eval(temp_ret));
% temp_func =@(xc) temp_ret_func(xc) - PHI3(xc)+PHI6(xc)+I3 + m * J_n * phi_c(xc, n_idx);
% temp_func(0.5)

main;
close all;
temp_omega = 2 * omega_m;
x = 1;
temp_target = PHI1(x)-PHI4(x)-I1;
eval(subs(-temp_omega^2 * psi_ck_sym - diff(psi_ck_sym, xc, 2) + lambda^2 * int(psi_ck_sym, xc, 0, 1) + lambda^2 * (H*tan(theta) - 1/2) * subs(psi_ck_sym, xc, 1), xc, x)) - temp_target
eval(subs(psi_ck_sym,xc,1) - subs(psi_bk_sym,xb,1) * cos(theta)^2)
eval(subs(psi_bk_sym_p,xb,0))
eval(subs(psi_bk_sym_p2,xb,1))
eval(subs(psi_bk_sym_p3,xb,1) - (k1 * subs(psi_ck_sym_p,xc,1) + k2 * int(psi_ck_sym, xc, 0, 1) + k3 * subs(psi_bk_sym, xb, 1) + k4 * subs(psi_bk_sym_p, xb, 1)))