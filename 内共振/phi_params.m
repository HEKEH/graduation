Ib = K * Ec * Ac^2 / Eb; %m4
beta1 = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);
k1 =  m / (beta1 * cos(theta));
k2 =  m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta));
k3 =  m * (- 2 * H * tan(theta) + 1)^2 * lambda^2 * cos(theta) / (4 * beta1);
k4 =  -m / (beta1 * cos(theta));

omegas_num = 8;
deter_func = @(omegas) determinant(omegas, beta1, theta, H, lambda, P, k1, k2, k3, k4);
omegas = get_omegas(deter_func, omegas_num);

alpha_b = sqrt(sqrt(P^2 + 4 * beta1 * omegas.^2) / 2 / beta1 + P / 2 / beta1);
beta_b = sqrt(sqrt(P^2 + 4 * beta1 * omegas.^2) / 2 / beta1 - P / 2 / beta1);
Gamma1 = lambda^2 * sin(omegas) ./ omegas ./ (omegas.^2 - lambda^2) + (cos(omegas) + lambda^2 * sin(omegas) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);
Gamma2 = lambda^2 * (1 - cos(omegas)) ./ omegas ./ (omegas.^2 - lambda^2) + (sin(omegas) + lambda^2 * (1 - cos(omegas)) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);

C = zeros(6, omegas_num);
C(1, :) = -1 * (Gamma2 * (cos(theta))^2 .* (alpha_b.^2 + beta_b.^2) .* (alpha_b .* cos(alpha_b) .* sinh(beta_b) - beta_b .* sin(alpha_b) .* cosh(beta_b))) ./ (Gamma1 .* sin(omegas) - Gamma2 .* (cos(omegas) - 1) + sin(omegas)) ./ (alpha_b.^2 .* cos(alpha_b)+ beta_b.^2 .* cosh(beta_b));
C(2, :) = 1 * ((Gamma1 + 1) * (cos(theta))^2 .* (alpha_b.^2 + beta_b.^2) .* (alpha_b .* cos(alpha_b) .* sinh(beta_b) - beta_b .* sin(alpha_b) .* cosh(beta_b))) ./ (Gamma1 .* sin(omegas) - Gamma2 .* (cos(omegas) - 1) + sin(omegas)) ./ (alpha_b.^2 .* cos(alpha_b)+ beta_b.^2 .* cosh(beta_b));
C(3, :) = 1 * alpha_b .* beta_b .* (alpha_b .* sin(alpha_b) + beta_b .* sinh(beta_b)) ./ (alpha_b.^2 .* cos(alpha_b)+ beta_b.^2 .* cosh(beta_b));
C(4, :) = -1 * beta_b;
C(5, :) = - 1 * alpha_b .* beta_b .* (alpha_b .* sin(alpha_b) + beta_b .* sinh(beta_b)) ./ (alpha_b.^2 .* cos(alpha_b)+ beta_b.^2 .* cosh(beta_b));
C(6, :) = 1 * alpha_b;

phi_c = @(xc, idx) C(1, idx)' .* (cos(omegas(idx)' * xc) + Gamma1(idx)') + C(2, idx)' .* (sin(omegas(idx)' * xc) + Gamma2(idx)');
phi_b = @(xb, idx) C(3, idx)' .* cos(alpha_b(idx)' * xb) + C(4, idx)' .* sin(alpha_b(idx)' * xb) + C(5, idx)' .* cosh(beta_b(idx)' * xb) + C(6, idx)' .* sinh(beta_b(idx)' * xb);

for i = 1: omegas_num
    C(:, i) = C(:, i) / sqrt(integral(@(x)m * phi_c(x, i) .^ 2 + cos(theta)^3 * phi_b(x, i) .^ 2, 0, 1));
end
phi_c = @(xc, idx) C(1, idx)' .* (cos(omegas(idx)' * xc) + Gamma1(idx)') + C(2, idx)' .* (sin(omegas(idx)' * xc) + Gamma2(idx)');
phi_b = @(xb, idx) C(3, idx)' .* cos(alpha_b(idx)' * xb) + C(4, idx)' .* sin(alpha_b(idx)' * xb) + C(5, idx)' .* cosh(beta_b(idx)' * xb) + C(6, idx)' .* sinh(beta_b(idx)' * xb);
phi_c_p = @(xc, idx) omegas(idx) * (-C(1, idx)' * sin(omegas(idx) * xc) +C(2, idx)' * cos(omegas(idx) * xc));
phi_c_p2 = @(xc, idx) omegas(idx)^2 * (-C(1, idx)' * cos(omegas(idx) * xc) - C(2, idx)' * sin(omegas(idx) * xc));
