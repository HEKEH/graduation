psi2_params_sym;

alpha_m2 = sqrt(sqrt(P.^2 + 16 .* beta1 .* omegas.^2) ./ 2 ./ beta1 + P ./ 2 ./ beta1);
beta_m2 = sqrt(sqrt(P.^2 + 16 .* beta1 .* omegas.^2) ./ 2 ./ beta1 - P ./ 2 ./ beta1);

Gamma_c21 = lambda.^2 .* (omegas .* (2.* H .* tan(theta) - 1) .* cos(2 .* omegas) + sin(2 .* omegas)) ./ (8 .* omegas .^ 3 - lambda.^2 .* omegas .* (2 .* H .* tan(theta)+1));
Gamma_c22 = lambda.^2 .* (omegas .* (2 .* H .* tan(theta) - 1) .* sin(2 .* omegas)+2 .* sin(omegas) .^ 2) ./ (8 .* omegas .^ 3 - lambda.^2  .* omegas .* (2 .* H .* tan(theta) + 1));
Ups_c21 = (H.*lambda.^2.*C(1, :).*(C(1, :).*(2.*sin(omegas)+omegas.*(-cos(omegas)+Gamma1+2.*H.*(cos(omegas)+Gamma1).*tan(theta)))+C(2, :).*(2-2.*cos(omegas)+omegas.*(-sin(omegas)+Gamma2+2.*H.*(sin(omegas)+Gamma2).*tan(theta)))))./(6.*omegas);
Ups_c22 = (H.*lambda.^2.*C(2, :).*(C(1, :).*(2.*sin(omegas)+omegas.*(-cos(omegas)+Gamma1+2.*H.*(cos(omegas)+Gamma1).*tan(theta)))+C(2, :).*(2-2.*cos(omegas)+omegas.*(-sin(omegas)+Gamma2+2.*H.*(sin(omegas)+Gamma2).*tan(theta)))))./(6.*omegas);

temp1 = (lambda.^2.*(2.*sin(omegas)+cos(omegas).*omegas.*(-1+2.*H.*tan(theta))))./(8.*omegas.^3-lambda.^2.*omegas.*(1+2.*H.*tan(theta)));
temp2 = (lambda.^2.*(2-2.*cos(omegas)+sin(omegas).*omegas.*(-1+2.*H.*tan(theta))))./(8.*omegas.^3-lambda.^2.*omegas.*(1+2.*H.*tan(theta)));
temp3 = (H.*lambda.^2.*omegas.*(4.*sin(omegas).^2.*C(1, :).*C(2, :)+C(1, :).^2.*((sin(2.*omegas)-2.*omegas))-C(2, :).^2.*((sin(2.*omegas)+2.*omegas))))./(2*(-8.*omegas.^2+lambda.^2.*(1+2.*H.*tan(theta))));
Ups_c23 = temp1 .* Ups_c21 + temp2 .* Ups_c22 + temp3;

C_c21 = zeros(1, omegas_num);
C_c22 = zeros(1, omegas_num);
C_b21 = zeros(1, omegas_num);
C_b22 = zeros(1, omegas_num);
C_b23 = zeros(1, omegas_num);
C_b24 = zeros(1, omegas_num);


for i = 1: omegas_num
    Ups_ck1 = Ups_c21(i);
    Ups_ck2 = Ups_c22(i);
    Ups_ck3 = Ups_c23(i);
    omega_m = omegas(i);
    alpha_bk = alpha_m2(i);
    beta_bk = beta_m2(i);
    Gamma_k1 = Gamma_c21(i);
    Gamma_k2 = Gamma_c22(i);
    temp = double(eval(C_k_sym));
    
    C_c21(i) = temp(1);
    C_c22(i) = temp(2);
    C_b21(i) = temp(3);
    C_b22(i) = temp(4);
    C_b23(i) = temp(5);
    C_b24(i) = temp(6);
end

psi_c2 = @(xc, idx) Ups_c21(idx) * cos(omegas(idx) * xc) + Ups_c22(idx) * sin(omegas(idx) * xc) + Ups_c23(idx) + C_c21(idx) * (cos(2 * omegas(idx) * xc) + Gamma_c21(idx)) + C_c22(idx) * (sin(2 * omegas(idx) * xc) + Gamma_c22(idx));
psi_b2 = @(xb, idx) C_b21(idx) * cos(alpha_m2(idx) * xb) + C_b22(idx) * sin(alpha_m2(idx) * xb) + C_b23(idx) * cosh(beta_m2(idx) * xb) + C_b24(idx) * sinh(beta_m2(idx) * xb);
psi_c2_p = @(xc, idx) omegas(idx) * (-Ups_c21(idx) * sin(omegas(idx) * xc) + Ups_c22(idx) * cos(omegas(idx) * xc)) + 2 * omegas(idx) * (-C_c21(idx) * (sin(2 * omegas(idx) * xc)) + C_c22(idx) * (cos(2 * omegas(idx) * xc)));
psi_c2_p2 = @(xc, idx) -omegas(idx)^2 * (Ups_c21(idx) * cos(omegas(idx) * xc) + Ups_c22(idx) * sin(omegas(idx) * xc)) - (2 * omegas(idx))^2 * (C_c21(idx) * (cos(2 * omegas(idx) * xc)) + C_c22(idx) * (sin(2 * omegas(idx) * xc)));