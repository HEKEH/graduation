psi3_params_sym;
alpha_m3 = repmat(sqrt(P / beta1), 1, omegas_num);

Ups_c31 = -((H.*lambda.^2.*C(1, :).*(C(1, :).*(2.*sin(omegas)+omegas.*(-cos(omegas)+Gamma1+2.*H.*(cos(omegas)+Gamma1).*tan(theta)))+C(2, :).*(2-2.*cos(omegas)+omegas.*(-sin(omegas)+Gamma2+2.*H.*(sin(omegas)+Gamma2).*tan(theta)))))./(2.*omegas));
Ups_c32 = -((H.*lambda.^2.*C(2, :).*(C(1, :).*(2.*sin(omegas)+omegas.*(-cos(omegas)+Gamma1+2.*H.*(cos(omegas)+Gamma1).*tan(theta)))+C(2, :).*(2-2.*cos(omegas)+omegas.*(-sin(omegas)+Gamma2+2.*H.*(sin(omegas)+Gamma2).*tan(theta)))))./(2.*omegas));

temp1 = -((2.*sin(omegas)+cos(omegas).*omegas.*(-1+2.*H.*tan(theta)))./(omegas.*(1+2.*H.*tan(theta))));
temp2 = (2*(-1+cos(omegas))-sin(omegas).*omegas.*(-1+2.*H.*tan(theta)))./(omegas.*(1+2.*H.*tan(theta)));
temp3 = (H.*omegas.*(4.*sin(omegas).^2.*C(1, :).*C(2, :)+C(1, :).^2.*((sin(2.*omegas)-2.*omegas))-C(2, :).^2.*((sin(2.*omegas)+2.*omegas))))./(2+4.*H.*tan(theta));
Ups_c33 = temp1 .* Ups_c31 + temp2 .* Ups_c32 + temp3;

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
psi_c3_p2 = @(xc, idx) omegas(idx)^2 * (-Ups_c31(idx) * cos(omegas(idx) * xc) - Ups_c32(idx) * sin(omegas(idx) * xc)) + 2 * C_c32(idx);