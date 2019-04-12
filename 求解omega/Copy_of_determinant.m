function ret = determinant(omega, m, beta1, theta, H, lambda, P)
    k1 =  - m / (beta1 * cos(theta)) * (1 / (4 * H^2) + 1);
    k2 =  m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta));
    k3 =  m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta)) * (H * sin(theta) * cos(theta) - (cos(theta)) ^ 2 / 2);
    k4 =  m * (2 * H + tan(theta)) / (2 * H * beta1 * cos(theta));

    alpha_b = sqrt(sqrt(P^2 + 4 * beta1 * omega^2) / 2 / beta1 + P / 2 / beta1);
    beta_b = sqrt(sqrt(P^2 + 4 * beta1 * omega^2) / 2 / beta1 - P / 2 / beta1);
    Gamma1 = lambda^2 * sin(omega) / omega / (omega^2 - lambda^2) + 1 / (2*omega^2-lambda^2) * (cos(omega) + lambda^2 * sin(omega) / omega / (omega^2 - lambda^2)) / (1 / lambda^2 / (2 * H * tan(theta) - 1) - 1 / 2 / omega^2 - lambda^4 / 2 / omega^2 / (omega^2 - lambda^2));
    Gamma2 = lambda^2 * (1 - cos(omega)) / omega / (omega^2 - lambda^2) + 1 / (2*omega^2-lambda^2) * (sin (omega) + lambda^2 * (1 - cos(omega)) / omega / (omega^2 - lambda^2)) / (1 / lambda^2 / (2 * H * tan(theta) - 1) - 1 / 2 / omega^2 - lambda^4 / 2 / omega^2 / (omega^2 - lambda^2));
    
    Lambda1 = (alpha_b^2 + beta_b^2) * (alpha_b * cos(alpha_b) * sinh(beta_b) - beta_b * sin(alpha_b) * cosh(beta_b));
    Lambda2 = k3 * Lambda1 + alpha_b * beta_b * (-alpha_b ^ 4 - beta_b^4 + alpha_b * beta_b * sin(alpha_b) * sinh(beta_b) * (-alpha_b^2 + beta_b^2-2 * k4) + cos(alpha_b) * cosh(beta_b) * (k4 * alpha_b^2 - beta_b^2 * (2 * alpha_b^2 + k4)) - k4 * alpha_b^2 + k4 * beta_b^2);
    ret = Lambda1 * (cos(theta))^2 / omega * (k1 * omega^2 *(-Gamma2 * sin(omega) + Gamma1 * cos(omega) + cos(omega)) + k2 * (Gamma2 * (sin(omega) - omega) + (Gamma1 + 1) * (cos(omega) - 1))) - 2 * Lambda2 * sin(omega / 2) * (Gamma2 * sin(omega / 2) + (Gamma1 + 1) * cos(omega / 2));
end
