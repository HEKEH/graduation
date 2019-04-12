psi2_params_sym;
Configs;
K = 100;
phi_params;

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

x = 0: 1e-2: 1;
yc = zeros(omegas_num, length(x));
yb = zeros(omegas_num, length(x));
for i = 1: omegas_num
    yc(i, :) = lc * psi_c2(x, i);
    yb(i, :) = lb * psi_b2(x, i);
end

for i = 1: omegas_num
    temp = 1 / sqrt(integral(@(x)psi_c2(x, i) .^ 2 , 0, 1) + integral(@(x)psi_b2(x, i) .^ 2, 0, 1));
    yci = yc(i, :) * temp / 20;
    ybi = yb(i, :) * temp / 20;
%     yci = yc(i, :) / 1000;
%     ybi = yb(i, :) / 1000;
    hold off
    subplot(4, 2, i)
    a = plot(lb * x + (ybi(end) * x * sin(theta) * cos(theta) - yci * sin(theta)), ...
        lb * x(end: -1: 1) * tan(theta) - (yci * cos(theta) + ybi(end) * x * sin(theta) * sin(theta)), 'red', 'linewidth', 0.75);
    hold on
    plot(lb * x, lb * x(end: -1: 1) * tan(theta), 'red--', 'linewidth', 0.75)
    hold on
    b = plot(lb * x, - ybi, 'blue', 'linewidth', 1.5);
    hold on
    plot(lb * x , zeros(1, length(x)), 'blue--', 'linewidth', 1.5);
    
    %axis off;
    axis equal;
    axis([0, lb * 1.1, -0.1 * lb, 1.1 * lb * tan(theta)])
    title(['{\it\Psi}{\it_c_2}、{\it\Psi}{\it_b_2}图({\itm}=', num2str(i), ')'], 'FontSize', 10);
    legend([a, b], {'{\it\Psi}{\it_c_2}', '{\it\Psi}{\it_b_2}'})
    set(gca, 'XTick', []); % 清除X轴的记号点
    set(gca, 'XGrid','off'); % X轴的网格
    set(gca, 'YTick', []); 
    set(gca, 'YGrid','off');
    set(gca, 'Position', [mod(i + 1, 2) * (0.5 - 1 / 60) + 1 / 30, floor((8 - i) / 2) * 0.25 + 0.018, 0.45, 0.195]);
end
set (gcf,'Position',[20,-200,600,1000], 'color','w') 

table = zeros(10, 8);
table(1, :) = Ups_c21;
table(2, :) = Ups_c22;
table(3, :) = Ups_c23;
table(4, :) = C_c21;
table(5, :) = C_c22;
table(6, :) = Gamma_c21;
table(7, :) = Gamma_c22;
table(8, :) = C_b21;
table(9, :) = C_b22;
table(10, :) = C_b23;
table(11, :) = C_b24;
table(12, :) = alpha_m2;
table(13, :) = beta_m2;
table = table';

lambda_ck = zeros(1, omegas_num);
lambda_bk = zeros(1, omegas_num);
for i = 1: omegas_num
    e_c = m * integral(@(x)psi_c2(x, i) .^ 2 , 0, 1);
    e_b = cos(theta)^3 * integral(@(x)psi_b2(x, i) .^ 2 , 0, 1);
    lambda_ck(i) = e_c / (e_c + e_b);
    lambda_bk(i) = e_b / (e_c + e_b);
end
lambda_bck = [lambda_ck; lambda_bk];