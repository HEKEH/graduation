mu_c = 1e-3;
mu_b = 1e-3;

psi_cks = cell(6, 1);
psi_bks = cell(6, 1);
psi_ck_ps = cell(6, 1);
psi_bk_ps = cell(6, 1);
psi_ck_p2s = cell(6, 1);
psi_bk_p2s = cell(6, 1);
psi_ck_ints = cell(6, 1);

for p_idx = 1: 6
    this_folder = ['../../resonance21_', num2str(p_idx)];
    run([this_folder, '/main.m']);
    psi_cks{p_idx} = psi_ck;
    psi_bks{p_idx} = psi_bk;
    psi_ck_ps{p_idx} = psi_ck_p;
    psi_bk_ps{p_idx} = psi_bk_p;
    psi_ck_p2s{p_idx} = psi_ck_p2;
    psi_bk_p2s{p_idx} = psi_bk_p2;
    psi_ck_ints{p_idx} = psi_ck_int;
end

close all 
%参数
Fc = @(xc) 1e-7 * ones(size(xc, 1), size(xc, 2));
Fb = @(xb) 1e-7 * m * ones(size(xb, 1), size(xb, 2));
sigma1 = omega_m - 2 * omega_n;
%参数

PHI_c31 = @(xc)-(1/2)*H*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{3}(1)-H*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{5}(1)+H*lambda^2*phi_c_p2(xc,m_idx)*integral(@(xc)psi_cks{3}(xc), 0, 1)+2*H*lambda^2*phi_c_p2(xc,m_idx)*integral(@(xc)psi_cks{5}(xc), 0, 1)-2*H*lambda^2*integral(@(xc)phi_c_p(xc,m_idx).*psi_ck_ps{3}(xc), 0, 1)-4*H*lambda^2*integral(@(xc)phi_c_p(xc,m_idx).*psi_ck_ps{5}(xc), 0, 1)-1/2*H*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{3}(xc)+H*lambda^2*psi_ck_p2s{3}(xc)*integral(@(xc)phi_c(xc,m_idx), 0, 1)-H*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{5}(xc)+2*H*lambda^2*psi_ck_p2s{5}(xc)*integral(@(xc)phi_c(xc,m_idx), 0, 1)+H^2*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{3}(1)*tan(theta)+2*H^2*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{5}(1)*tan(theta)+H^2*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{3}(xc)*tan(theta)+2*H^2*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{5}(xc)*tan(theta);

PHI_c32 = @(xc)-(1/2)*H*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{2}(1)-1/2*H*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{4}(1)-H*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{6}(1)+H*lambda^2*phi_c_p2(xc,n_idx)*integral(@(xc)psi_cks{2}(xc), 0, 1)+H*lambda^2*phi_c_p2(xc,n_idx)*integral(@(xc)psi_cks{4}(xc), 0, 1)+2*H*lambda^2*phi_c_p2(xc,m_idx)*integral(@(xc)psi_cks{6}(xc), 0, 1)-2*H*lambda^2*integral(@(xc)phi_c_p(xc,n_idx).*psi_ck_ps{2}(xc), 0, 1)-2*H*lambda^2*integral(@(xc)phi_c_p(xc,n_idx).*psi_ck_ps{4}(xc), 0, 1)-4*H*lambda^2*integral(@(xc)phi_c_p(xc,m_idx).*psi_ck_ps{6}(xc), 0, 1)-1/2*H*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{2}(xc)+H*lambda^2*psi_ck_p2s{2}(xc)*integral(@(xc)phi_c(xc,n_idx), 0, 1)-1/2*H*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{4}(xc)+H*lambda^2*psi_ck_p2s{4}(xc)*integral(@(xc)phi_c(xc,n_idx), 0, 1)-H*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{6}(xc)+2*H*lambda^2*psi_ck_p2s{6}(xc)*integral(@(xc)phi_c(xc,m_idx), 0, 1)+H^2*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{2}(1)*tan(theta)+H^2*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{4}(1)*tan(theta)+2*H^2*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{6}(1)*tan(theta)+H^2*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{2}(xc)*tan(theta)+H^2*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{4}(xc)*tan(theta)+2*H^2*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{6}(xc)*tan(theta);

PHI_c33 = @(xc)-(1/2)*H*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{2}(1)-1/2*H*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{4}(1)-H*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{5}(1)+H*lambda^2*phi_c_p2(xc,m_idx)*integral(@(xc)psi_cks{2}(xc), 0, 1)+H*lambda^2*phi_c_p2(xc,m_idx)*integral(@(xc)psi_cks{4}(xc), 0, 1)+2*H*lambda^2*phi_c_p2(xc,n_idx)*integral(@(xc)psi_cks{5}(xc), 0, 1)-2*H*lambda^2*integral(@(xc)phi_c_p(xc,m_idx).*psi_ck_ps{2}(xc), 0, 1)-2*H*lambda^2*integral(@(xc)phi_c_p(xc,m_idx).*psi_ck_ps{4}(xc), 0, 1)-4*H*lambda^2*integral(@(xc)phi_c_p(xc,n_idx).*psi_ck_ps{5}(xc), 0, 1)-1/2*H*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{2}(xc)+H*lambda^2*psi_ck_p2s{2}(xc)*integral(@(xc)phi_c(xc,m_idx), 0, 1)-1/2*H*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{4}(xc)+H*lambda^2*psi_ck_p2s{4}(xc)*integral(@(xc)phi_c(xc,m_idx), 0, 1)-H*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{5}(xc)+2*H*lambda^2*psi_ck_p2s{5}(xc)*integral(@(xc)phi_c(xc,n_idx), 0, 1)+H^2*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{2}(1)*tan(theta)+H^2*lambda^2*phi_c_p2(xc,m_idx)*psi_cks{4}(1)*tan(theta)+2*H^2*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{5}(1)*tan(theta)+H^2*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{2}(xc)*tan(theta)+H^2*lambda^2*phi_c(1,m_idx)*psi_ck_p2s{4}(xc)*tan(theta)+2*H^2*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{5}(xc)*tan(theta);

PHI_c34 = @(xc)-(1/2)*H*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{1}(1)-H*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{6}(1)+H*lambda^2*phi_c_p2(xc,n_idx)*integral(@(xc)psi_cks{1}(xc), 0, 1)+2*H*lambda^2*phi_c_p2(xc,n_idx)*integral(@(xc)psi_cks{6}(xc), 0, 1)-2*H*lambda^2*integral(@(xc)phi_c_p(xc,n_idx).*psi_ck_ps{1}(xc), 0, 1)-4*H*lambda^2*integral(@(xc)phi_c_p(xc,n_idx).*psi_ck_ps{6}(xc), 0, 1)-1/2*H*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{1}(xc)+H*lambda^2*psi_ck_p2s{1}(xc)*integral(@(xc)phi_c(xc,n_idx), 0, 1)-H*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{6}(xc)+2*H*lambda^2*psi_ck_p2s{6}(xc)*integral(@(xc)phi_c(xc,n_idx), 0, 1)+H^2*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{1}(1)*tan(theta)+2*H^2*lambda^2*phi_c_p2(xc,n_idx)*psi_cks{6}(1)*tan(theta)+H^2*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{1}(xc)*tan(theta)+2*H^2*lambda^2*phi_c(1,n_idx)*psi_ck_p2s{6}(xc)*tan(theta);

Gamma_cm1 = @(xc)PHI_c31(xc);
Gamma_cm2 = @(xc)PHI_c32(xc) - (2 * psi_cks{1}(xc) * omega_m - phi_c(xc, m_idx) * m * J_m / 2 / omega_m) * m * J_n / omega_n;
Gamma_cm3 = @(xc)-psi_cks{1}(xc) * omega_m + phi_c(xc, m_idx) * m * J_m / 2 / omega_m;
Gamma_bm1 = @(xb)0;
Gamma_bm2 = @(xb)- (2 * psi_bks{1}(xb) * omega_m - phi_b(xb, m_idx) * m * J_m / 2 / omega_m) * m * J_n / omega_n;
Gamma_bm3 = @(xb)-psi_bks{1}(xb) * omega_m + phi_b(xb, m_idx) * m * J_m / 2 / omega_m;

Gamma_cn1 = @(xc)PHI_c33(xc) + m * J_n * psi_cks{2}(xc) - phi_c(xc, n_idx) * (m * J_n / 2 / omega_n)^2;
Gamma_cn2 = @(xc)PHI_c34(xc) - m * J_m / omega_m * (psi_cks{2}(xc) * omega_n - phi_c(xc, n_idx) * m * J_n / 4 / omega_n);
Gamma_cn3 = @(xc)psi_cks{2}(xc) * omega_n - phi_c(xc, n_idx) * m * J_n / 2 / omega_n;
Gamma_bn1 = @(xb)m * J_n * psi_bks{2}(xb) - phi_b(xb, n_idx) * (m * J_n / 2 / omega_n)^2;
Gamma_bn2 = @(xb)-m * J_m / omega_m * (psi_bks{2}(xb) * omega_n - phi_b(xb, n_idx) * m * J_n / 4 / omega_n);
Gamma_bn3 = @(xb)psi_bks{2}(xb) * omega_n - phi_b(xb, n_idx) * m * J_n / 2 / omega_n;

pm = zeros(1, 4);
pn = zeros(1, 6);
pm(1) = m * integral(@(xc)phi_c(xc, m_idx) .^ 2, 0, 1);
pm(2) = m * integral(@(xc)phi_c(xc, m_idx) .* Gamma_cm1(xc), 0, 1) + cos(theta)^3 * integral(@(xb)phi_b(xb, m_idx) .* Gamma_bm1(xb), 0, 1);
pm(3) = m * integral(@(xc)phi_c(xc, m_idx) .* Gamma_cm2(xc), 0, 1) + cos(theta)^3 * integral(@(xb)phi_b(xb, m_idx) .* Gamma_bm2(xb), 0, 1);
pm(4) = m * integral(@(xc)phi_c(xc, m_idx) .* Gamma_cm3(xc), 0, 1) + cos(theta)^3 * integral(@(xb)phi_b(xb, m_idx) .* Gamma_bm3(xb), 0, 1);
% pm(5) = 1/2 * m * integral(@(xc)phi_c(xc, m_idx) .* Fc(xc), 0, 1);
% pm(6) = 1/2 * cos(theta)^3 * integral(@(xb)phi_b(xb, m_idx) .* Fb(xb), 0, 1);

pn(1) = m * integral(@(xc)phi_c(xc, n_idx) .^ 2, 0, 1);
pn(2) = m * integral(@(xc)phi_c(xc, n_idx) .* Gamma_cn1(xc), 0, 1) + cos(theta)^3 * integral(@(xb)phi_b(xb, n_idx) .* Gamma_bn1(xb), 0, 1);
pn(3) = m * integral(@(xc)phi_c(xc, n_idx) .* Gamma_cn2(xc), 0, 1) + cos(theta)^3 * integral(@(xb)phi_b(xb, n_idx) .* Gamma_bn2(xb), 0, 1);
pn(4) = m * integral(@(xc)phi_c(xc, n_idx) .* Gamma_cn3(xc), 0, 1) + cos(theta)^3 * integral(@(xb)phi_b(xb, n_idx) .* Gamma_bn3(xb), 0, 1);
pn(5) = 1/2 * m * integral(@(xc)phi_c(xc, n_idx) .* Fc(xc), 0, 1);
pn(6) = 1/2 * cos(theta)^3 * integral(@(xb)phi_b(xb, n_idx) .* Fb(xb), 0, 1);
