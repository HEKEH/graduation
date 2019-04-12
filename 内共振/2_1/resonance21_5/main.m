params;
C_k = double(eval(C_k_sym));
psi_ck_sym = eval(psi_ck_known_sym) + C_k' * eval(psi_ck_uncertain_sym');
psi_bk_sym = eval(psi_bk_known_sym) + C_k' * eval(psi_bk_uncertain_sym');
psi_ck_sym_p =  eval(psi_ck_known_sym_p) + C_k' * eval(psi_ck_uncertain_sym_p');
psi_bk_sym_p = eval(psi_bk_known_sym_p) + C_k' * eval(psi_bk_uncertain_sym_p');
psi_ck_sym_p2 =  eval(psi_ck_known_sym_p2) + C_k' * eval(psi_ck_uncertain_sym_p2');
psi_bk_sym_p2 = eval(psi_bk_known_sym_p2) + C_k' * eval(psi_bk_uncertain_sym_p2');
psi_bk_sym_p3 = eval(psi_bk_known_sym_p3) + C_k' * eval(psi_bk_uncertain_sym_p3');

psi_ck = matlabFunction(psi_ck_sym);
psi_bk = matlabFunction(psi_bk_sym);
psi_ck_p = matlabFunction(psi_ck_sym_p);
psi_bk_p = matlabFunction(psi_bk_sym_p);
psi_ck_p2 = matlabFunction(psi_ck_sym_p2);
psi_bk_p2 = matlabFunction(psi_bk_sym_p2);
psi_bk_p3 = matlabFunction(psi_bk_sym_p3);
psi_ck_int = double(eval(psi_ck_known_sym_int)) + C_k' * double(eval(psi_ck_uncertain_sym_int))';
e_c = m * double(int(psi_ck_sym^2, xc, 0, 1));
e_b = cos(theta)^3 * double(int(psi_bk_sym^2, xb, 0, 1));


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
title('{\it\Psi}{\it_c_5}、{\it\Psi}{\it_b_5}图', 'FontSize', 10);
%title(['{\it\Psi}{\it_c_m}、{\it\Psi}{\it_b_m}图 (放大',num2str(temp / 20),'倍)'] ,'FontSize', 10);
legend([a, b], {'{\it\Psi}{\it_c_5}', '{\it\Psi}{\it_b_5}'})
set(gca, 'XTick', []); % 清除X轴的记号点
set(gca, 'XGrid','off'); % X轴的网格
set(gca, 'YTick', []); 
set(gca, 'YGrid','off');
set(gca, 'Position', [0.05, 0.01, 0.9, 0.94]);
set(gcf,'Position',[20,20,800,550], 'color','w') 
% 
table = zeros(1, 15);
table(1) = omega_m;
table(2) = omega_n;
table(3) = Ups_ck1;
table(4) = Ups_ck2;
table(5) = Ups_ck3;
table(6) = nan;
table(7) = nan;
table(8) = alpha_bk;
table(9) = nan;
table(10: 15) = C_k;
table(1: 2) = [];

lambda_ck = e_c / (e_c + e_b);
lambda_bk = e_b / (e_c + e_b);
lambda_bck = [lambda_ck; lambda_bk];