alpha_m3 = repmat(sqrt(P / beta1), 1, omegas_num);

Ups_c31 = -(C(1,:) .* H .* lambda.^2) ./ omegas .* (C(1,:) .* (omegas .* (Gamma1 + H .* tan(theta) .* (Gamma1+cos(omegas))) + sin(omegas)) + C(2,:) .* (omegas .* (Gamma2 + H .* tan(theta) .* (Gamma2 + sin(omegas)))-cos(omegas)+1));
Ups_c32 = -(C(2,:) .* H .* lambda.^2) ./ omegas .* (C(1,:) .* (omegas .* (Gamma1 + H .* tan(theta) .* (Gamma1+cos(omegas))) + sin(omegas)) + C(2,:) .* (omegas .* (Gamma2 + H .* tan(theta) .* (Gamma2 + sin(omegas)))-cos(omegas)+1));

Ups_c33 = -H .* omegas .* (2 .* (C(1, :) .^2 + C(2, :) .^2) .* omegas - (4 .* C(1, :) .* C(2, :) .* sin(omegas).^2 + (C(1, :) .^2 - C(2, :) .^2) .* sin(2 .* omegas))) ./ 2 ./ (2 .* H .* tan(theta) + 1)...
    + ((1- 2 .* H .* tan(theta)) .* (Ups_c31 .* cos(omegas) + Ups_c32 .* sin(omegas)) - H) ./ (2 .* H .* tan(theta) + 1)...
    - 2 .* (Ups_c31 .* sin(omegas) + Ups_c32 .* (1 - cos(omegas))) ./ omegas ./ (2 .* H .* tan(theta) + 1);

C_c31 = zeros(1, omegas_num);
C_c32 = zeros(1, omegas_num);
C_b31 = zeros(1, omegas_num);
C_b32 = zeros(1, omegas_num);
C_b33 = zeros(1, omegas_num);
C_b34 = zeros(1, omegas_num);

AA = - 2 * H * tan(theta) / (2 * H * tan(theta) + 1);
BB = - (lambda^2 * (6 * H * tan(theta) - 1) - 12) / (3 * lambda^2 * (2 * H * tan(theta) + 1));
Gamma_k1 = AA;
Gamma_k2 = BB;

for i = 1: omegas_num
    Ups_ck1 = Ups_c31(i);
    Ups_ck2 = Ups_c32(i);
    Ups_ck3 = Ups_c33(i);
    omega_m = omegas(i);
    alpha_bk = alpha_m3(i);
    temp = double(eval(C_k_sym));
    
    C_c31(i) = temp(1);
    C_c32(i) = temp(2);
    C_b31(i) = temp(3);
    C_b32(i) = temp(4);
    C_b33(i) = temp(5);
    C_b34(i) = temp(6);
end

psi_c3 = @(xc, idx) Ups_c31(idx) * cos(omegas(idx) * xc) + Ups_c32(idx) * sin(omegas(idx) * xc) +Ups_c33(idx) + C_c31(idx) * (xc + AA) + C_c32(idx) * (xc.^2 + BB);
psi_b3 = @(xb, idx) C_b31(idx) * cos(alpha_m3(idx) * xb) + C_b32(idx) * sin(alpha_m3(idx) * xb) + C_b33(idx) * xb + C_b34(idx);
psi_c3_p = @(xc, idx) omegas(idx) * (-Ups_c31(idx) * sin(omegas(idx) * xc) + Ups_c32(idx) * cos(omegas(idx) * xc)) + C_c31(idx) + 2 * C_c32(idx) * xc;
psi_c3_p2 = @(xc, idx) omegas(idx)^2 * (-Ups_c31(idx) * cos(omegas(idx) * xc) - Ups_c32(idx) * sin(omegas(idx) * xc)) + C_c31(idx) + 2 * C_c32(idx);