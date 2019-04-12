matr_A = sym(zeros(6, 6));
vec_b = sym(zeros(6, 1));

matr_A(1, :) = subs(psi_ck_uncertain_sym, xc, 0);
vec_b(1) = -subs(psi_ck_known_sym, xb, 0);

matr_A(2, :) = subs(psi_bk_uncertain_sym, xb, 0);
vec_b(2) = -subs(psi_bk_known_sym, xb, 0);

matr_A(3, :) = subs(psi_bk_uncertain_sym_p, xb, 0);
vec_b(3) = -subs(psi_bk_known_sym_p, xb, 0);

matr_A(4, :) = subs(psi_bk_uncertain_sym_p2, xb, 1); 
vec_b(4) = -subs(psi_bk_known_sym_p2, xb, 1);

matr_A(5, :) = subs(psi_ck_uncertain_sym, xc, 1) - subs(psi_bk_uncertain_sym, xb, 1) * cos(theta)^2;
vec_b(5) = -subs(psi_ck_known_sym, xc, 1) + subs(psi_bk_known_sym, xb, 1) * cos(theta)^2;

matr_A(6, :) = -subs(psi_bk_uncertain_sym_p3, xb, 1) + k1 * subs(psi_ck_uncertain_sym_p, xc, 1) + k2 * psi_ck_uncertain_sym_int + k3 * subs(psi_bk_uncertain_sym, xb, 1) + k4 * subs(psi_bk_uncertain_sym_p, xb, 1);
vec_b(6) = subs(psi_bk_known_sym_p3, xb, 1) - k1 * subs(psi_ck_known_sym_p, xc, 1) - k2 * subs(psi_ck_known_sym_int, xc, 1) - k3 * subs(psi_bk_known_sym, xb, 1) - k4 * subs(psi_bk_known_sym_p, xb, 1);

