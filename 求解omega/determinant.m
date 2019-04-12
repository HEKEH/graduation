function ret = determinant(omegas, beta1, theta, H, lambda, P, k1, k2, k3, k4)
    alpha_b = sqrt(sqrt(P^2 + 4 * beta1 * omegas.^2) / 2 / beta1 + P / 2 / beta1);
    beta_b = sqrt(sqrt(P^2 + 4 * beta1 * omegas.^2) / 2 / beta1 - P / 2 / beta1);
    Gamma1 = lambda^2 * sin(omegas) ./ omegas ./ (omegas.^2 - lambda^2) + (cos(omegas) + lambda^2 * sin(omegas) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);
    Gamma2 = lambda^2 * (1 - cos(omegas)) ./ omegas ./ (omegas.^2 - lambda^2) + (sin(omegas) + lambda^2 * (1 - cos(omegas)) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);
%     Lambda1 = (alpha_b.^2 + beta_b.^2) .* (alpha_b .* cos(alpha_b) .* sinh(beta_b) - beta_b .* sin(alpha_b) .* cosh(beta_b));
%     Lambda2 = k3 * Lambda1 + alpha_b .* beta_b .* (-alpha_b .^ 4 - beta_b.^4 + alpha_b .* beta_b .* sin(alpha_b) .* sinh(beta_b) .* (-alpha_b.^2 + beta_b.^2-2 * k4) + cos(alpha_b) .* cosh(beta_b) .* (k4 * alpha_b.^2 - beta_b.^2 .* (2 * alpha_b.^2 + k4)) - k4 * alpha_b.^2 + k4 .* beta_b.^2);
    Lambda1 = (alpha_b.^2 + beta_b.^2) .* (alpha_b .* cos(alpha_b) .* sinh(beta_b) - beta_b .* sin(alpha_b) .* cosh(beta_b));
    Lambda2 = k3 .* Lambda1 + alpha_b .* beta_b .* (-alpha_b .^ 4 - beta_b.^4 + alpha_b .* beta_b .* sin(alpha_b) .* sinh(beta_b) .* (-alpha_b.^2 + beta_b.^2-2 .* k4) + cos(alpha_b) .* cosh(beta_b) .* (k4 .* alpha_b.^2 - beta_b.^2 .* (2 .* alpha_b.^2 + k4)) - k4 .* alpha_b.^2 + k4 .* beta_b.^2);
    ret = Lambda1 * (cos(theta))^2 ./ omegas .* (-k1 * omegas.^2 .*(Gamma2 .* sin(omegas) + Gamma1 .* cos(omegas) + cos(omegas)) + k2 * (Gamma2 .* (sin(omegas) - omegas) + (Gamma1 + 1) .* (cos(omegas) - 1))) - 2 * Lambda2 .* sin(omegas ./ 2) .* (Gamma2 .* sin(omegas / 2) + (Gamma1 + 1) .* cos(omegas / 2));
end
