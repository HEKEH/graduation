params;
psi_ck_known=matlabFunction(eval(psi_ck_known_sym));
psi_ck_known_p=matlabFunction(eval(psi_ck_known_sym_p));
psi_ck_known_p_2=matlabFunction(eval(psi_ck_known_sym_p2));
psi_ck_known_int=double(eval(psi_ck_known_sym_int));

psi_bk_known=matlabFunction(eval(psi_bk_known_sym));
psi_bk_known_p=matlabFunction(eval(psi_bk_known_sym_p));
psi_bk_known_p2=matlabFunction(eval(psi_bk_known_sym_p2));
psi_bk_known_p3=matlabFunction(eval(psi_bk_known_sym_p3));

psi_ck_uncertain=matlabFunction(eval(psi_ck_uncertain_sym));
psi_ck_uncertain_p=matlabFunction(eval(psi_ck_uncertain_sym_p));
psi_ck_uncertain_p2=matlabFunction(eval(psi_ck_uncertain_sym_p2));
psi_ck_uncertain_int=double(eval(psi_ck_uncertain_sym_int));
psi_bk_uncertain=matlabFunction(eval(psi_bk_uncertain_sym));
psi_bk_uncertain_p=matlabFunction(eval(psi_bk_uncertain_sym_p));
psi_bk_uncertain_p2=matlabFunction(eval(psi_bk_uncertain_sym_p2));
psi_bk_uncertain_p3=matlabFunction(eval(psi_bk_uncertain_sym_p3));

matr_A = zeros(6, 6);
vec_b = zeros(6, 1);

matr_A(1, :) = psi_ck_uncertain(0);
vec_b(1) = -psi_ck_known(0);

matr_A(2, :) = psi_bk_uncertain(0);
vec_b(2) = -psi_bk_known(0);

matr_A(3, :) = psi_bk_uncertain_p(0);
vec_b(3) = -psi_bk_known_p(0);

matr_A(4, :) = psi_bk_uncertain_p2(1); 
vec_b(4) = -psi_bk_known_p2(1);

matr_A(5, :) = psi_ck_uncertain(1) - psi_bk_uncertain(1) * cos(theta)^2;
vec_b(5) = -psi_ck_known(1) + psi_bk_known(1) * cos(theta)^2;

matr_A(6, :) = -psi_bk_uncertain_p3(1) + k1 * psi_ck_uncertain_p(1) + k2 * psi_ck_uncertain_int + k3 * psi_bk_uncertain(1) + k4 * psi_bk_uncertain_p(1);
vec_b(6) = psi_bk_known_p3(1) - k1 * psi_ck_known_p(1) - k2 * psi_ck_known_int(1) - k3 * psi_bk_known(1) - k4 * psi_bk_known_p(1);

C_k = matr_A \ vec_b;

psi_ck = @(xc) psi_ck_known(xc) + C_k' * psi_ck_uncertain(xc')';
psi_bk = @(xb) psi_bk_known(xb) + C_k' * psi_bk_uncertain(xb')';

dx = 1e-2;
x = 0: dx: 1;
len = length(x);
yc = zeros(1, len);
yb = zeros(1, len);
for i = 1: len
    yc(i) = psi_ck(x(i));
    yb(i) = psi_bk(x(i));
end

temp = 1 / sqrt(sum(yc.^2) * dx + sum(yb.^2) * dx);
yc = lc * yc * temp / 20;
yb = lb * yb * temp / 20;
%     yc = yc / 1000;
%     yb = yb / 1000;
hold off
a = plot(lb * x +(yb(end) * x * sin(theta) * cos(theta) - yc * sin(theta)),... 
    lb * x(end: -1: 1) * tan(theta) -(yc * cos(theta) + yb(end) * x * sin(theta) * sin(theta)), 'red', 'linewidth', 0.75);
hold on
plot(lb * x, lb * x(end: -1: 1) * tan(theta), 'red--', 'linewidth', 0.75)
hold on
b = plot(lb * x, - yb, 'blue', 'linewidth', 1.5);
hold on
plot(lb * x , zeros(1, length(x)), 'blue--', 'linewidth', 1.5);

%axis off;
axis equal;
axis([0, lb * 1.1, -0.1 * lb, 1.1 * lb * tan(theta)])
title('{\it\Psi}{\it_c_1}、{\it\Psi}{\it_b_1}图', 'FontSize', 10);
%title(['{\it\Psi}{\it_c_m}、{\it\Psi}{\it_b_m}图 (放大',num2str(temp / 20),'倍)'] ,'FontSize', 10);
legend([a, b], {'{\it\Psi}{\it_c_1}', '{\it\Psi}{\it_b_1}'})
set(gca, 'XTick', []); % 清除X轴的记号点
set(gca, 'XGrid','off'); % X轴的网格
set(gca, 'YTick', []); 
set(gca, 'YGrid','off');
set(gca, 'Position', [0.05, 0.01, 0.9, 0.94]);
set(gcf,'Position',[20,20,800,550], 'color','w') 
% 
table = zeros(1, 19);
table(1) = omega_m;
table(2) = omega_n;
table(3) = Ups_ck1;
table(4) = Ups_ck2;
table(5) = Ups_ck3;
table(6) = Ups_ck4;
table(7) = Ups_ck5;
table(8) = alpha_bk;
table(9) = beta_bk;
table(10) = Ups_bk1;
table(11) = Ups_bk2;
table(12) = Ups_bk3;
table(13) = Ups_bk4;
table(14: 19) = C_k;

e_c = m * double(int((eval(psi_ck_known_sym) + C_k' * eval(psi_ck_uncertain_sym'))^2, xc, 0, 1));
e_b = cos(theta)^3 * double(int((eval(psi_bk_known_sym) + C_k' * eval(psi_bk_uncertain_sym'))^2, xb, 0, 1));
lambda_ck = e_c /(e_c + e_b);
lambda_bk = e_b /(e_c + e_b);
lambda_bck = [lambda_ck; lambda_bk];