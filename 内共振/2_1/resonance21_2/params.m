func_sym;
run('../../Configs.m');
run('../../Kmn_param.m');
run('../../phi_params.m');
C_idx = n_idx;%要改
Gamma_k1 = Gamma1(C_idx);%要改
Gamma_k2 = Gamma2(C_idx);%要改

Gamma1_m = Gamma1(m_idx);
Gamma2_m = Gamma2(m_idx);
Gamma1_n = Gamma1(n_idx);
Gamma2_n = Gamma2(n_idx);
C_m = C(:, m_idx);
C_n = C(:, n_idx);
omega_m = omegas(m_idx);
omega_n = omegas(n_idx);


%--------------------------要改-----------------------------
alpha_bk = sqrt(sqrt(P^2 + 4 * beta1 * omega_n^2) / 2 / beta1 + P / 2 / beta1);%要改
beta_bk = sqrt(sqrt(P^2 + 4 * beta1 * omega_n^2) / 2 / beta1 - P / 2 / beta1);%要改
%--------------------------要改-----------------------------

PHI1 = @(xc)H*lambda^2*phi_c_p2(xc, m_idx)*integral(@(xc)phi_c(xc, m_idx),0, 1);
PHI2 = @(xc)H*lambda^2*phi_c_p2(xc, n_idx)*integral(@(xc)phi_c(xc, n_idx),0, 1);
PHI3 = @(xc)H*lambda^2*(phi_c_p2(xc, m_idx)*integral(@(xc)phi_c(xc, n_idx),0, 1)+phi_c_p2(xc, n_idx)*integral(@(xc)phi_c(xc, m_idx),0, 1));
PHI4 = @(xc)H*lambda^2*(1/2-H*tan(theta))*phi_c_p2(xc, m_idx)*phi_c(1, m_idx);
PHI5 = @(xc)H*lambda^2*(1/2-H*tan(theta))*phi_c_p2(xc, n_idx)*phi_c(1, n_idx);
PHI6 = @(xc)H*lambda^2*(1/2-H*tan(theta))*(phi_c_p2(xc, m_idx)*phi_c(1, n_idx)+phi_c_p2(xc, n_idx)*phi_c(1, m_idx));
I1 = H*lambda^2*integral(@(xc)phi_c_p(xc, m_idx).^2,0, 1);
I2 = H*lambda^2*integral(@(xc)phi_c_p(xc, n_idx).^2,0, 1);
I3 = 2*H*lambda^2*integral(@(xc)phi_c_p(xc, m_idx).*phi_c_p(xc, n_idx),0, 1);

J_m = integral(@(xc)phi_c(xc, m_idx) .* (PHI2(xc) - PHI5(xc) - I2), 0, 1);
J_n = integral(@(xc)phi_c(xc, n_idx) .* (PHI3(xc) - PHI6(xc) - I3), 0, 1);

%--------------------------要改-----------------------------
Ups_ck1 = -((C_n(2)*(2*m*J_n*omega_m+H*lambda^2*omega_n^2*(C_m(1)*(2*sin(omega_m)+omega_m*(cos(omega_m)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma1_m))+C_m(2)*(2-2*cos(omega_m)+omega_m*(sin(omega_m)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma2_m)))))/(4*omega_m*omega_n));

Ups_ck2 = (C_n(1)*(2*m*J_n*omega_m+H*lambda^2*omega_n^2*(C_m(1)*(2*sin(omega_m)+omega_m*(cos(omega_m)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma1_m))+C_m(2)*(2-2*cos(omega_m)+omega_m*(sin(omega_m)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma2_m)))))/(4*omega_m*omega_n);

Ups_ck3 = (H*lambda^2*omega_m^2*C_m(1)*(C_n(1)*(2*sin(omega_n)+omega_n*(cos(omega_n)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma1_n))+C_n(2)*(2-2*cos(omega_n)+omega_n*(sin(omega_n)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma2_n))))/(2*(-omega_m^2*omega_n+omega_n^3));

Ups_ck4 = (H*lambda^2*omega_m^2*C_m(2)*(C_n(1)*(2*sin(omega_n)+omega_n*(cos(omega_n)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma1_n))+C_n(2)*(2-2*cos(omega_n)+omega_n*(sin(omega_n)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma2_n))))/(2*(-omega_m^2*omega_n+omega_n^3));

temp1 = (lambda^2*(2-2*sin(omega_n)*omega_n+cos(omega_n)*(-2+omega_n^2*((1-2*H*tan(theta))))))/(-2*omega_n^4+lambda^2*omega_n^2*(1+2*H*tan(theta)));

temp2 = (lambda^2*(2*cos(omega_n)*omega_n+sin(omega_n)*(-2+omega_n^2*((1-2*H*tan(theta))))))/(-2*omega_n^4+lambda^2*omega_n^2*(1+2*H*tan(theta)));

temp3 = -((lambda^2*(2*sin(omega_m)+cos(omega_m)*omega_m*(-1+2*H*tan(theta))))/(omega_m*(-2*omega_n^2+lambda^2*(1+2*H*tan(theta)))));

temp4 = (lambda^2*(2*(-1+cos(omega_m))-sin(omega_m)*omega_m*(-1+2*H*tan(theta))))/(omega_m*(-2*omega_n^2+lambda^2*(1+2*H*tan(theta))));

temp5 = (1/((omega_m^2-omega_n^2)*(-2*omega_n^2+lambda^2*(1+2*H*tan(theta)))))*2*(-2*H*lambda^2*omega_m*omega_n^2*(C_n(1)*C_m(2)+sin(omega_m)*C_m(1)*(cos(omega_n)*C_n(1)+sin(omega_n)*C_n(2))-cos(omega_m)*C_m(2)*(cos(omega_n)*C_n(1)+sin(omega_n)*C_n(2)))+m*J_n*omega_n^2*(C_n(1)*Gamma1_n+C_n(2)*Gamma2_n)+omega_m^2*((2*H*lambda^2*omega_n*((sin(omega_n)*C_n(1)*((cos(omega_m)*C_m(1)+sin(omega_m)*C_m(2)))+((C_m(1)-cos(omega_n)*((cos(omega_m)*C_m(1)+sin(omega_m)*C_m(2)))))*C_n(2)))-m*J_n*((C_n(1)*Gamma1_n+C_n(2)*Gamma2_n)))));

Ups_ck5 = temp1 * Ups_ck1 + temp2 * Ups_ck2 + temp3 * Ups_ck3 + temp4 * Ups_ck4 + temp5;

Ups_bk1 = C_n(4)/(2 *beta1* alpha_bk*(alpha_bk^2 + beta_bk^2)) * (-m * J_n);
Ups_bk2 = -(C_n(3)/(2 *beta1* alpha_bk*(alpha_bk^2 + beta_bk^2))) * (-m * J_n);
Ups_bk3 = C_n(6)/(2 *beta1* beta_bk*(alpha_bk^2 + beta_bk^2)) * (-m * J_n);
Ups_bk4 = C_n(5)/(2 *beta1* beta_bk*(alpha_bk^2 + beta_bk^2)) * (-m * J_n);
%--------------------------要改-----------------------------