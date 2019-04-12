function ret = determinant2(omegas, beta1, theta, H, lambda, P, k1, k2, k3, k4)
alpha_b = sqrt(sqrt(P^2 + 4 * beta1 * omegas.^2) / 2 / beta1 + P / 2 / beta1);
    beta_b = sqrt(sqrt(P^2 + 4 * beta1 * omegas.^2) / 2 / beta1 - P / 2 / beta1);
    Gamma1 = lambda^2 * sin(omegas) ./ omegas ./ (omegas.^2 - lambda^2) + (cos(omegas) + lambda^2 * sin(omegas) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);
    Gamma2 = lambda^2 * (1 - cos(omegas)) ./ omegas ./ (omegas.^2 - lambda^2) + (sin(omegas) + lambda^2 * (1 - cos(omegas)) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);
%     Lambda1 = (alpha_b.^2 + beta_b.^2) .* (alpha_b .* cos(alpha_b) .* sinh(beta_b) - beta_b .* sin(alpha_b) .* cosh(beta_b));
%     Lambda2 = k3 * Lambda1 + alpha_b .* beta_b .* (-alpha_b .^ 4 - beta_b.^4 + alpha_b .* beta_b .* sin(alpha_b) .* sinh(beta_b) .* (-alpha_b.^2 + beta_b.^2-2 * k4) + cos(alpha_b) .* cosh(beta_b) .* (k4 * alpha_b.^2 - beta_b.^2 .* (2 * alpha_b.^2 + k4)) - k4 * alpha_b.^2 + k4 .* beta_b.^2);
%     ret = Lambda1 * (cos(theta))^2 ./ omegas .* (-k1 * omegas.^2 .*(Gamma2 .* sin(omegas) + Gamma1 .* cos(omegas) + cos(omegas)) + k2 * (Gamma2 .* (sin(omegas) - omegas) + (Gamma1 + 1) .* (cos(omegas) - 1))) - 2 * Lambda2 .* sin(omegas ./ 2) .* (Gamma2 .* sin(omegas / 2) + (Gamma1 + 1) .* cos(omegas / 2));
    len = length(omegas);
    ret = zeros(1,len);
    for i = 1: len
        ret(i) = det([1 + Gamma1(i), Gamma2(i), 0, 0, 0, 0;
                0, 0, 1, 0, 1, 0;
                0, 0, 0, alpha_b(i), 0, beta_b(i);
                0, 0, -alpha_b(i)^2 * cos(alpha_b(i)), -alpha_b(i)^2 * sin(alpha_b(i)), beta_b(i)^2 * cosh(beta_b(i)), beta_b(i)^2 * sinh(beta_b(i));
                cos(omegas(i)) + Gamma1(i), sin(omegas(i)) + Gamma2(i), -cos(alpha_b(i)) * cos(theta)^2, -sin(alpha_b(i)) * cos(theta)^2, -cosh(beta_b(i)) * cos(theta)^2, -sinh(beta_b(i)) * cos(theta)^2;
                -k1 * omegas(i) * sin(omegas(i)) + k2 * (sin(omegas(i)) / omegas(i) + Gamma1(i)),...
                k1 * omegas(i) * cos(omegas(i)) + k2 * ((1 - cos(omegas(i))) / omegas(i) + Gamma2(i)),...
                -alpha_b(i)^3 * sin(alpha_b(i)) + k3 * cos(alpha_b(i)) - k4 * alpha_b(i) * sin(alpha_b(i)),...
                alpha_b(i)^3 * cos(alpha_b(i)) + k3 * sin(alpha_b(i)) + k4 * alpha_b(i) * cos(alpha_b(i)),...
                -beta_b(i)^3 * sinh(beta_b(i)) + k3 * cosh(beta_b(i)) + k4 * beta_b(i) * sinh(beta_b(i)),...
                -beta_b(i)^3 * cosh(beta_b(i)) + k3 * sinh(beta_b(i)) + k4 * beta_b(i) * cosh(beta_b(i));
            ]);
    end
end
