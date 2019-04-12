Configs;
K = 100;
Ib = K * Ec * Ac^2 / Eb; %m4
beta1 = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);
k1 =  m / (beta1 * cos(theta));
k2 =  m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta));
k3 =  m * (- 2 * H * tan(theta) + 1)^2 * lambda^2 * cos(theta) / (4 * beta1);
k4 =  -m / (beta1 * cos(theta));

omegas_num = 8;
omegas = get_omegas(K, lb, Eb, lc, Ec, Ac, Hc, mb, mc, m, theta, H, lambda, P, 20);
omegas(omegas_num + 1: length(omegas)) = [];

alpha_b = sqrt(sqrt(P^2 + 4 * beta1 * omegas.^2) / 2 / beta1 + P / 2 / beta1);
beta_b = sqrt(sqrt(P^2 + 4 * beta1 * omegas.^2) / 2 / beta1 - P / 2 / beta1);
Gamma1 = lambda^2 * sin(omegas) ./ omegas ./ (omegas.^2 - lambda^2) + (cos(omegas) + lambda^2 * sin(omegas) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);
Gamma2 = lambda^2 * (1 - cos(omegas)) ./ omegas ./ (omegas.^2 - lambda^2) + (sin(omegas) + lambda^2 * (1 - cos(omegas)) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);

C = zeros(6, omegas_num);
% for i = 1: omegas_num
%     temp_mtx = [1 + Gamma1(i), Gamma2(i), 0, 0, 0, 0;
%                 0, 0, 1, 0, 1, 0;
%                 0, 0, 0, alpha_b(i), 0, beta_b(i);
%                 0, 0, -alpha_b(i)^2 * cos(alpha_b(i)), -alpha_b(i)^2 * sin(alpha_b(i)), beta_b(i)^2 * cosh(beta_b(i)), beta_b(i)^2 * sinh(beta_b(i));
%                 cos(omegas(i)) + Gamma1(i), sin(omegas(i)) + Gamma2(i), -cos(alpha_b(i)) * cos(theta)^2, -sin(alpha_b(i)) * cos(theta)^2, -cosh(beta_b(i)) * cos(theta)^2, -sinh(beta_b(i)) * cos(theta)^2;
%                 ];
%     C(:, i) = null(temp_mtx, 'r');
% end
C(1, :) = -1 * (Gamma2 * (cos(theta))^2 .* (alpha_b.^2 + beta_b.^2) .* (alpha_b .* cos(alpha_b) .* sinh(beta_b) - beta_b .* sin(alpha_b) .* cosh(beta_b))) ./ (Gamma1 .* sin(omegas) - Gamma2 .* (cos(omegas) - 1) + sin(omegas)) ./ (alpha_b.^2 .* cos(alpha_b)+ beta_b.^2 .* cosh(beta_b));
C(2, :) = 1 * ((Gamma1 + 1) * (cos(theta))^2 .* (alpha_b.^2 + beta_b.^2) .* (alpha_b .* cos(alpha_b) .* sinh(beta_b) - beta_b .* sin(alpha_b) .* cosh(beta_b))) ./ (Gamma1 .* sin(omegas) - Gamma2 .* (cos(omegas) - 1) + sin(omegas)) ./ (alpha_b.^2 .* cos(alpha_b)+ beta_b.^2 .* cosh(beta_b));
C(3, :) = 1 * alpha_b .* beta_b .* (alpha_b .* sin(alpha_b) + beta_b .* sinh(beta_b)) ./ (alpha_b.^2 .* cos(alpha_b)+ beta_b.^2 .* cosh(beta_b));
C(4, :) = -1 * beta_b;
C(5, :) = - 1 * alpha_b .* beta_b .* (alpha_b .* sin(alpha_b) + beta_b .* sinh(beta_b)) ./ (alpha_b.^2 .* cos(alpha_b)+ beta_b.^2 .* cosh(beta_b));
C(6, :) = 1 * alpha_b;

phi_c0 = @(xc, idx) C(1, idx)' .* (cos(omegas(idx)' * xc) + Gamma1(idx)') + C(2, idx)' .* (sin(omegas(idx)' * xc) + Gamma2(idx)');
phi_b0 = @(xb, idx) C(3, idx)' .* cos(alpha_b(idx)' * xb) + C(4, idx)' .* sin(alpha_b(idx)' * xb) + C(5, idx)' .* cosh(beta_b(idx)' * xb) + C(6, idx)' .* sinh(beta_b(idx)' * xb);

CC = zeros(omegas_num, 1);
for i = 1: omegas_num
    CC(i) = 1 / sqrt(m * integral(@(x)phi_c0(x, i) .^ 2 , 0, 1) + cos(theta)^3 * integral(@(x)phi_b0(x, i) .^ 2, 0, 1));
end
C = C .* repmat(CC', 6, 1);
x = 0: 1e-2: 1;
phi_c = @(xc, idx) CC(idx) .* phi_c0(xc, idx);
phi_b = @(xb, idx) CC(idx) .* phi_b0(xb, idx);
yc = lc .* phi_c(x, 1: omegas_num);
yb = lb .* phi_b(x, 1: omegas_num);

for i = 1: omegas_num
    temp = 1 / sqrt(integral(@(x)phi_c(x, i) .^ 2 , 0, 1) + integral(@(x)phi_b(x, i) .^ 2, 0, 1));
    yci = yc(i, :) * temp / 20;
    ybi = yb(i, :) * temp / 20;
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
    title(['第', num2str(i), '阶振型(\omega=', num2str(round(omegas(i), 4)), ')'], 'FontSize', 11);
    legend([a, b], {'索', '梁'})
    set(gca, 'XTick', []); % 清除X轴的记号点
    set(gca, 'XGrid','off'); % X轴的网格
    set(gca, 'YTick', []); 
    set(gca, 'YGrid','off');
    set(gca, 'Position', [mod(i + 1, 2) * (0.5 - 1 / 60) + 1 / 30, floor((8 - i) / 2) * 0.25 + 0.023, 0.45, 0.2]);
end
set (gcf,'Position',[20,20,600,1000], 'color','w') 

%局部化因子
lambda_ck = zeros(1, omegas_num);
lambda_bk = zeros(1, omegas_num);
for i = 1: omegas_num
    e_c = m * integral(@(x)phi_c0(x, i) .^ 2 , 0, 1);
    e_b = cos(theta)^3 * integral(@(x)phi_b0(x, i) .^ 2 , 0, 1);
    lambda_ck(i) = e_c / (e_c + e_b);
    lambda_bk(i) = e_b / (e_c + e_b);
end
lambda_bck = [lambda_ck; lambda_bk];
table=[omegas; alpha_b; beta_b; Gamma1; Gamma2; CC']';
table = [table, C'];