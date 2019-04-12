alpha_m2 = sqrt(sqrt(P.^2 + 16 .* beta1 .* omegas.^2) ./ 2 ./ beta1 + P ./ 2 ./ beta1);
beta_m2 = sqrt(sqrt(P.^2 + 16 .* beta1 .* omegas.^2) ./ 2 ./ beta1 - P ./ 2 ./ beta1);

Gamma_c21 = lambda.^2 .* (omegas .* (2.* H .* tan(theta) - 1) .* cos(2 .* omegas) + sin(2 .* omegas)) ./ (8 .* omegas .^ 3 - lambda.^2 .* omegas .* (2 .* H .* tan(theta)+1));
Gamma_c22 = lambda.^2 .* (omegas .* (2 .* H .* tan(theta) - 1) .* sin(2 .* omegas)+2 .* sin(omegas) .^ 2) ./ (8 .* omegas .^ 3 - lambda.^2  .* omegas .* (2 .* H .* tan(theta) + 1));
Ups_c21 = (C(1,:) .* H .* lambda.^2) ./ (3 .* omegas) .* (C(1,:) .* (omegas .* (Gamma1 + H .* tan(theta) .* (Gamma1+cos(omegas))) + sin(omegas)) + C(2,:) .* (omegas .* (Gamma2+H .* tan(theta) .* (Gamma2 + sin(omegas)))-cos(omegas)+1));
Ups_c22 = ((C(2,:) .* H .* lambda.^2) ./(3 .* omegas)) .* (C(1,:) .* (omegas .* (Gamma1 + H .* tan(theta) .* (Gamma1+cos(omegas))) + sin(omegas)) + C(2,:) .* (omegas .* (Gamma2+H .* tan(theta) .* (Gamma2 + sin(omegas)))-cos(omegas)+1));

Ups_c23=(4 .* lambda.^2 .* (Ups_c22 + Ups_c21 .* sin(omegas) - Ups_c22 .* cos(omegas))) ./ (16 .* omegas .^ 3 - 2 .* lambda.^2  .* omegas .* (2 .* H .* tan(theta) + 1))...
    -(H .* lambda.^2  .* omegas .* (4 .* C(1,:) .* C(2,:) .* sin(omegas) .^ 2 - C(2,:) .^ 2 .* (2 .* omegas+sin(2 .* omegas)) + C(1,:).^2 .* (sin(2 .* omegas)-2 .* omegas))) ./ (16 .* omegas .^ 2- 2 .* lambda.^2 .* (2 .* H .* tan(theta)+1))...
    +(lambda.^2) .* ((2 * H * tan(theta) - 1) .* (Ups_c22 .* sin(omegas)+Ups_c21 .* cos(omegas)) + H) ./ (8 .* omegas.^2-lambda.^2 .* (2 .* H .* tan(theta)+1));

C_c21 = zeros(1, omegas_num);
C_c22 = zeros(1, omegas_num);
C_b21 = zeros(1, omegas_num);
C_b22 = zeros(1, omegas_num);
C_b23 = zeros(1, omegas_num);
C_b24 = zeros(1, omegas_num);

matr_A = zeros(6, 6);
vec_b = zeros(6, 1);
for i = 1: omegas_num
    matr_A(1, 1) = (1 + Gamma_c21(i));
    matr_A(1, 2) = Gamma_c22(i);
    vec_b(1) = -Ups_c21(i)-Ups_c23(i);
    
    matr_A(2, 3) = 1;
    matr_A(2, 5) = 1;
    
    matr_A(3, 4) = alpha_m2(i);
    matr_A(3, 6) = beta_m2(i);    
    
    matr_A(4, 3) = -alpha_m2(i)^2 * cos(alpha_m2(i)); 
    matr_A(4, 4) = -alpha_m2(i)^2 * sin(alpha_m2(i));
    matr_A(4, 5) = +beta_m2(i)^2 * cosh(beta_m2(i));
    matr_A(4, 6) = +beta_m2(i)^2 * sinh(beta_m2(i));
    
    matr_A(5, 1) = (Gamma_c21(i) + cos(2 * omegas(i)));
    matr_A(5, 2) = (Gamma_c22(i) + sin(2 * omegas(i)));
    matr_A(5, 3) = - cos(theta)^2 * cos(alpha_m2(i)); 
    matr_A(5, 4) = - cos(theta)^2 * sin(alpha_m2(i));
    matr_A(5, 5) = - cos(theta)^2 * cosh(beta_m2(i));
    matr_A(5, 6) = - cos(theta)^2 * sinh(beta_m2(i));
    vec_b(5) = -Ups_c21(i) * cos(omegas(i)) - Ups_c22(i) * sin (omegas(i)) - Ups_c23(i);
    
    matr_A(6, 1) = (k2 * Gamma_c21(i) + k1 * (Gamma_c21(i) + cos(2 * omegas(i))) + (k2 * sin(omegas(i)) * cos(omegas(i)))/omegas(i));
    matr_A(6, 2) = (k2 * Gamma_c22(i) + k1 * (Gamma_c22(i) + sin(2 * omegas(i))) + k2 * sin(omegas(i))^2/omegas(i));
    matr_A(6, 3) = (k4 * alpha_m2(i) * (-sin(alpha_m2(i))) + k3 * cos(alpha_m2(i)) + (-alpha_m2(i)^3) * sin(alpha_m2(i)));
    matr_A(6, 4) = (k3 * sin(alpha_m2(i)) + k4 * alpha_m2(i) * cos(alpha_m2(i)) + alpha_m2(i)^3 * cos(alpha_m2(i)));
    matr_A(6, 5) = (k4 * beta_m2(i) * sinh(beta_m2(i)) + k3 * cosh(beta_m2(i)) + (-beta_m2(i)^3) * sinh(beta_m2(i)));
    matr_A(6, 6) = (k3 * sinh(beta_m2(i)) + k4 * beta_m2(i) * cosh(beta_m2(i)) + (-beta_m2(i)^3) * cosh(beta_m2(i)));
    vec_b(6) = -(k2 * sin(omegas(i))/omegas(i) + k1 * cos(omegas(i))) * Ups_c21(i) + (-(k2/omegas(i)) - k1 * sin(omegas(i)) + k2 * cos(omegas(i))/omegas(i)) * Ups_c22(i) - (k1 + k2) * Ups_c23(i);

    temp = matr_A \ vec_b;
    
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