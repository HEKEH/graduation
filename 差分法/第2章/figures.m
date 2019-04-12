run('./Configs.m');
run('./MySettings.m');
folders_num = 3;%3
folders_idx = 1: folders_num;
t_all = cell(1, folders_num);
wc_mid = cell(1, folders_num);
wc_end = cell(1, folders_num);
wb_mid = cell(1, folders_num);
wb_end = cell(1, folders_num);
for i = folders_idx %组
    kth = i;
    condition = ['实例', num2str(kth)];
    folder = [condition, '\数据'];
    interval = 1;
    rows = zeros(5, 1);
    cols = zeros(5, 1);
    temp = csvread([folder, '\t',num2str(1),'.csv']); 
    [rows(1), cols(1)] = size(temp(1: interval: end, :));
    temp = csvread([folder, '\Wc',num2str(1),'.csv']); 
    [rows(2), cols(2)] = size(temp(1: interval: end, :));
    temp = csvread([folder, '\Vc',num2str(1),'.csv']); 
    [rows(3), cols(3)] = size(temp(1: interval: end, :));
    temp = csvread([folder, '\Wb',num2str(1),'.csv']); 
    [rows(4), cols(4)] = size(temp(1: interval: end, :));
    temp = csvread([folder, '\Vb',num2str(1),'.csv']); 
    [rows(5), cols(5)] = size(temp(1: interval: end, :));

    t = zeros(rows(1) * rounds, cols(1));
    Wc = zeros(rows(2) * rounds, cols(2));
    Vc = zeros(rows(3) * rounds, cols(3));
    Wb = zeros(rows(4) * rounds, cols(4));
    Vb = zeros(rows(5) * rounds, cols(5));

    for r = 1: rounds
        temp1 = csvread([folder, '\t',num2str(r),'.csv']); 
        temp2 = csvread([folder, '\Wc',num2str(r),'.csv']); 
        %temp3 = csvread([folder, '\Vc',num2str(r),'.csv']); 
        temp4 = csvread([folder, '\Wb',num2str(r),'.csv']); 
        %temp5 = csvread([folder, '\Vb',num2str(r),'.csv']); 

        t(rows(1) * (r - 1) + 1: rows(1) * r, :) = temp1(1: interval: end, :); 
        Wc(rows(2) * (r - 1) + 1: rows(2) * r, :) = temp2(1: interval: end, :); 
        %Vc(rows(3) * (r - 1) + 1: rows(3) * r, :) = temp3(1: interval: end, :); 
        Wb(rows(4) * (r - 1) + 1: rows(4) * r, :) = temp4(1: interval: end, :); 
        %Vb(rows(5) * (r - 1) + 1: rows(5) * r, :) = temp5(1: interval: end, :);
    end
    t_all{kth} = t;
    wc_mid{kth} = Wc(:, floor(Nc / 2 + 1));
    wc_end{kth} = Wc(:, Nc + 1);
    wb_mid{kth} = Wb(:, floor(Nb / 2 + 1));
    wb_end{kth} = Wb(:, Nb + 1);
end
%save figures.mat
w_all = {wc_mid, wc_end, wb_mid, wb_end};
points = {'索中点', '索端点', '梁中点', '梁端点'};
lines = {'-','-','-','-'};
fig_num = 1;
for i = folders_idx
    pL = zeros(1, 4);
    %lgd = cell(1, 4);
    figure(fig_num);
    fig_num = fig_num + 1;
    for j = 1: 4
        w_temp = w_all{j};
        pL(j) = plot(t_all{i}, w_temp{i},lines{j});
        hold on
        %lgd{i} = ['工况', num2str(3 * (i - 1) + k)];
    end
    hold off
    title(['实例', num2str(i),'振动时间历程曲线'],'FontSize', 14);%, 'FontSize', 10
    legend(pL, points,'Location','NorthWest');
%     if j == 1 || j == 2
%         w_name = 'w_c';
%     else
%         w_name = 'w_b';
%     end
    xlabel('无量纲时间t')
    ylabel('无量纲位移w_c/w_b')
    set (gcf,'Position',[50,50,1200,500], 'color','w') 
    set(gca, 'Position', [0.06, 0.09, 0.9, 0.85]);
end

lines = {'-','-.','--',':'};
start_times = [0, 200];
end_times = [25, 225];

for idx = 1: length(start_times)
    start_time = start_times(idx) + dt * dNt * interval;
    end_time = end_times(idx) + dt * dNt * interval;
    time_idxs = (start_time/ (dt * dNt * interval)): (end_time/ (dt * dNt * interval));
    for i = 1: folders_num
        pL = zeros(1, 4);
        %lgd = cell(1, 4);
        figure(fig_num);
        fig_num = fig_num + 1;
        for j = 1: 4
            w_temp = w_all{j};
            pL(j) = plot(t_all{i}(time_idxs), w_temp{i}(time_idxs), lines{j});
            hold on
            %lgd{i} = ['工况', num2str(3 * (i - 1) + k)];
        end
        hold off
        title(['实例', num2str(i),'振动时间历程曲线(片段', num2str(idx), ')'],'FontSize', 14);%, 'FontSize', 10
        legend(pL, points,'Location','NorthWest');
    %     if j == 1 || j == 2
    %         w_name = 'w_c';
    %     else
    %         w_name = 'w_b';
    %     end
        xlabel('无量纲时间t')
        ylabel('无量纲位移w_c/w_b')
        xlim([start_time, end_time-5])
        set (gcf,'Position',[50,50,1200,500], 'color','w') 
        set(gca, 'Position', [0.06, 0.09, 0.9, 0.85]);
    end
end

enlarge_times = 50;%放大倍数
Wc_real = Wc * lc;
Wb_real = Wb * lb;

start_time = 100.5;
end_time = 102.5;
pn = 8;
dtime = (end_time - start_time) / pn;
alltime = start_time: dtime: end_time - dtime;
time_idxs = round(alltime/ (dt * dNt * interval) + 1);

figure(fig_num);
fig_num = fig_num + 1;
for j = 1: pn
    i = time_idxs(j);
    subplot(4,2,j);
    a = plot(lb * (0: dxc: 1) + enlarge_times * (Wb_real(i, Nb) * (0: dxc: 1) * sin(theta) * cos(theta) - Wc_real(i, :) * sin(theta)), ...
        lb * (1: -dxc: 0) * tan(theta) - enlarge_times * (Wc_real(i, :) * cos(theta) + Wb_real(i, Nb) * (0: dxc: 1) * sin(theta) * sin(theta)),...
        'linewidth', 1.5,'color', 'r');
    hold on
    b = plot(lb * (0: dxb: 1), - enlarge_times * Wb_real(i, :), 'linewidth', 3,'color', 'b');
    axis equal
    axis([0, ceil(lb / 10) * 10, -10, Inf])
    title(['t: ', num2str(t(i))],'FontSize', 12);
    %xlim([0, ceil(lb / 10) * 10]);
    %legend([a, b], {'{\it\Psi}{\it_c_2}', '{\it\Psi}{\it_b_2}'})
    set(gca, 'XTick', []); % 清除X轴的记号点
    set(gca, 'XGrid','off'); % X轴的网格
    set(gca, 'YTick', []); 
    set(gca, 'YGrid','off');
    set(gca, 'Position', [mod(j + 1, 2) * (0.5 - 1 / 60) + 1 / 30, floor((8 - j) / 2) * 0.25 + 0.018, 0.45, 0.195]);
end
hold off
set (gcf,'Position',[100,-200,600,1000], 'color','w') 
