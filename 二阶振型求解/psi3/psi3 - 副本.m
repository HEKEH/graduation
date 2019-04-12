Configs;
K = 100;
phi_params;

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

matr_A = zeros(6, 6);
vec_b = zeros(6, 1);
AA = - 2 * H * tan(theta) / (2 * H * tan(theta) + 1);
BB = - (lambda^2 * (6 * H * tan(theta) - 1) - 12) / (3 * lambda^2 * (2 * H * tan(theta) + 1));


for i = 1: omegas_num
    matr_A(1, 1) = AA;
    matr_A(1, 2) = BB;
    vec_b(1) = -Ups_c31(i)-Ups_c33(i);
    
    matr_A(2, 3) = 1;
    matr_A(2, 6) = 1;
    
    matr_A(3, 4) = alpha_m3(i);
    matr_A(3, 5) = 1;    
    
    matr_A(4, 3) = -cos(alpha_m3(i)) * alpha_m3(i)^2;
    matr_A(4, 4) = -sin(alpha_m3(i)) * alpha_m3(i)^2;
    
    matr_A(5, 1) = 1 + AA;
    matr_A(5, 2) = 1 + BB;
    matr_A(5, 3) = -cos(theta)^2 * cos(alpha_m3(i)); 
    matr_A(5, 4) = -cos(theta)^2 * sin(alpha_m3(i));
    matr_A(5, 5) = -cos(theta)^2;
    matr_A(5, 6) = -cos(theta)^2;
    vec_b(5) = -Ups_c31(i) * cos(omegas(i)) - Ups_c32(i) * sin (omegas(i)) - Ups_c33(i);
  
    matr_A(6, 1) = k1 + k2 * (1/2 + AA);
    matr_A(6, 2) = 2 * k1 + k2 * (1/3 + BB);
    matr_A(6, 3) = -k4 * alpha_m3(i) * sin(alpha_m3(i)) + k3 * cos(alpha_m3(i)) - alpha_m3(i)^3 * sin(alpha_m3(i)); 
    matr_A(6, 4) = k4 * alpha_m3(i) * cos(alpha_m3(i)) + k3 * sin(alpha_m3(i)) + alpha_m3(i)^3 * cos(alpha_m3(i));
    matr_A(6, 5) = k3 + k4;
    matr_A(6, 6) = k3;
    vec_b(6) = -Ups_c31(i) * (k1 * cos(omegas(i)) + k2 * sin(omegas(i)) / omegas(i))- Ups_c32(i) * (k2 / omegas(i) + k1 * sin(omegas(i)) - k2 * cos(omegas(i)) / omegas(i))- (k1 + k2) * Ups_c33(i);

    temp = matr_A \ vec_b;
    
    C_c31(i) = temp(1);
    C_c32(i) = temp(2);
    C_b31(i) = temp(3);
    C_b32(i) = temp(4);
    C_b33(i) = temp(5);
    C_b34(i) = temp(6);
end
m_idx = 1;

psi_c3 = @(xc, idx) Ups_c31(idx) * cos(omegas(idx) * xc) + Ups_c32(idx) * sin(omegas(idx) * xc) +Ups_c33(idx) + C_c31(idx) * (xc + AA) + C_c32(idx) * (xc.^2 + BB);
psi_b3 = @(xb, idx) C_b31(idx) * cos(alpha_m3(idx) * xb) + C_b32(idx) * sin(alpha_m3(idx) * xb) + C_b33(idx) * xb + C_b34(idx);

x = 0: 1e-2: 1;
yc = zeros(omegas_num, length(x));
yb = zeros(omegas_num, length(x));
for i = 1: omegas_num
    yc(i, :) = lc * psi_c3(x, i);
    yb(i, :) = lb * psi_b3(x, i);
end

for i = 1: omegas_num
    temp = 1 / sqrt(integral(@(x)psi_c3(x, i) .^ 2 , 0, 1) + integral(@(x)psi_b3(x, i) .^ 2, 0, 1));
    yci = yc(i, :) * temp / 20;
    ybi = yb(i, :) * temp / 20;
%     yci = yc(i, :) / 10000;
%     ybi = yb(i, :) / 10000;
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
    title(['{\it\Psi}{\it_c_3}、{\it\Psi}{\it_b_3}图({\itm}=', num2str(i), ')'], 'FontSize', 10);
    legend([a, b], {'{\it\Psi}{\it_c_3}', '{\it\Psi}{\it_b_3}'})
    set(gca, 'XTick', []); % 清除X轴的记号点
    set(gca, 'XGrid','off'); % X轴的网格
    set(gca, 'YTick', []); 
    set(gca, 'YGrid','off');
    set(gca, 'Position', [mod(i + 1, 2) * (0.5 - 1 / 60) + 1 / 30, floor((8 - i) / 2) * 0.25 + 0.018, 0.45, 0.195]);
end
set (gcf,'Position',[20,20,600,1000], 'color','w') 

table = zeros(10, 8);
table(1, :) = Ups_c31;
table(2, :) = Ups_c32;
table(3, :) = Ups_c33;
table(4, :) = C_c31;
table(5, :) = C_c32;
table(6, :) = C_b31;
table(7, :) = C_b32;
table(8, :) = C_b33;
table(9, :) = C_b34;
table(10, :) = alpha_m3;
table = table';

lambda_ck = zeros(1, omegas_num);
lambda_bk = zeros(1, omegas_num);
for i = 1: omegas_num
    e_c = m * integral(@(x)psi_c3(x, i) .^ 2 , 0, 1);
    e_b = cos(theta)^3 * integral(@(x)psi_b3(x, i) .^ 2 , 0, 1);
    lambda_ck(i) = e_c / (e_c + e_b);
    lambda_bk(i) = e_b / (e_c + e_b);
end
lambda_bck = [lambda_ck; lambda_bk];