clear;
psi_c_3to6 = cell(6, 1);
psi_b_3to6 = cell(6, 1);
lambda_bcks = zeros(2, 6);

for p_idx = 3: 6
    this_folder = ['resonance21_', num2str(p_idx)];
    run([this_folder, '/main.m']);
    psi_c_3to6{p_idx} = psi_ck;
    psi_b_3to6{p_idx} = psi_bk;
    lambda_bcks(:, p_idx) = lambda_bck;
end
%clearvars -except psi_c_3to6 psi_b_3to6
close all 
dx = 1e-2;
x = 0: dx: 1;
len = length(x);
for p_idx = 3: 6
    yc = zeros(1, len);
    yb = zeros(1, len);
    for i = 1: len
        yc(i) = psi_c_3to6{p_idx}(x(i));
        yb(i) = psi_b_3to6{p_idx}(x(i));
    end

    temp = 1 / sqrt(sum(yc.^2) * dx + sum(yb.^2) * dx);
    yc = lc * yc * temp / 20;
    yb = lb * yb * temp / 20;
    %     yc = yc / 1000;
    %     yb = yb / 1000;
    subplot(2, 2, p_idx - 2);
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
    title(['{\it\Psi}{\it_c_',num2str(p_idx),'}、{\it\Psi}{\it_b_',num2str(p_idx),'}图'], 'FontSize', 14);
    %title(['{\it\Psi}{\it_c_m}、{\it\Psi}{\it_b_m}图 (放大',num2str(temp / 20),'倍)'] ,'FontSize', 10);
    l1 = legend([a, b], {['{\it\Psi}{\it_c_',num2str(p_idx),'}'], ['{\it\Psi}{\it_b_',num2str(p_idx),'}']});
    set(l1, 'FontSize',12)
    set(gca, 'XTick', []); % 清除X轴的记号点
    set(gca, 'XGrid','off'); % X轴的网格
    set(gca, 'YTick', []); 
    set(gca, 'YGrid','off');
    set(gca, 'Position', [mod(p_idx + 1, 2) * (0.5 - 1 / 60) + 1 / 30, floor((6 - p_idx) / 2) * 0.5 + 0.018, 0.45, 0.45]);
end
set (gcf,'Position',[20,20,1000,800], 'color','w') 