%parpool('local', 4)
clear;
if ~exist('data.mat', 'file')
    Configs;
    MySettings;
    folder = '数据';
    interval = 20;

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

    t_real = zeros(rows(1) * rounds, cols(1));
    Wc_real = zeros(rows(2) * rounds, cols(2));
    Vc_real = zeros(rows(3) * rounds, cols(3));
    Wb_real = zeros(rows(4) * rounds, cols(4));
    Vb_real = zeros(rows(5) * rounds, cols(5));
    
    for r = 1: rounds
        temp1 = csvread([folder, '\t',num2str(r),'.csv']); 
        temp2 = csvread([folder, '\Wc',num2str(r),'.csv']); 
        %temp3 = csvread([folder, '\Vc',num2str(r),'.csv']); 
        temp4 = csvread([folder, '\Wb',num2str(r),'.csv']); 
        %temp5 = csvread([folder, '\Vb',num2str(r),'.csv']); 

        t_real(rows(1) * (r - 1) + 1: rows(1) * r, :) = temp1(1: interval: end, :); 
        Wc_real(rows(2) * (r - 1) + 1: rows(2) * r, :) = temp2(1: interval: end, :); 
        %Vc_real(rows(3) * (r - 1) + 1: rows(3) * r, :) = temp3(1: interval: end, :); 
        Wb_real(rows(4) * (r - 1) + 1: rows(4) * r, :) = temp4(1: interval: end, :); 
        %Vb_real(rows(5) * (r - 1) + 1: rows(5) * r, :) = temp5(1: interval: end, :);
    end
    save data.mat
else
    if ~exist('Ac', 'var')
        load('data.mat')
    end
end

enlarge_times = 5e1;%放大倍数
for i = 1: 10: length(t_real)
    hold off
    plot(lb * (0: dxc: 1) + enlarge_times * (Wb_real(i, Nb) * (0: dxc: 1) * sin(theta) * cos(theta) - Wc_real(i, :) * sin(theta)), ...
        lb * (1: -dxc: 0) * tan(theta) - enlarge_times * (Wc_real(i, :) * cos(theta) + Wb_real(i, Nb) * (0: dxc: 1) * sin(theta) * sin(theta)),...
        'linewidth', 1.5);
    hold on
    plot(lb * (0: dxb: 1), - enlarge_times * Wb_real(i, :), 'linewidth', 3);
    axis equal
    axis([0, ceil(lb / 10) * 10, -10, Inf])
    title(['实际时间: ', num2str(t_real(i)), 's 无量纲时间: ', num2str(t_real(i) / lc * sqrt(Hc / mc)), ' (振动放大', num2str(enlarge_times), '倍)']);
    %xlim([0, ceil(lb / 10) * 10]);
    drawnow;
end

% hc = Ec * Ac / Le * ((sin(theta) - mc * g * lc * cos(theta)^2 / 2 / Hc) * Wb_real(:, Nb + 1) + mc * g * cos(theta) / Hc * Wc_real * [1/2; ones(Nc - 1, 1); 1/2] * dxc * lc);
% tan_phi = mc * g * lc * cos(theta) / 2 / Hc;
% Wblbp = (3 * Wb_real(:, Nb + 1) - 4 * Wb_real(:, Nb) + Wb_real(:, Nb - 1)) / (2 * dxb * lb);
% Wclcp = (3 * Wc_real(:, Nc + 1) - 4 * Wc_real(:, Nc) + Wc_real(:, Nc - 1)) / (2 * dxc * lc);
% Wblbp3 = (3 * Wb_real(:, Nb - 3) - 14 * Wb_real(:, Nb - 2) + 24 * Wb_real(:, Nb - 1) - 18 * Wb_real(:, Nb) + 5 * Wb_real(:, Nb + 1)) / (2 * (dxb * lb)^3);
% Wblbp3_theory = (hc * (sin(theta) - cos(theta) * tan_phi) - Hc * Wblbp * cos(theta) + Hc * Wclcp * cos(theta)) / Eb / Ib;
% Wblbp3 ./ Wblbp3_theory
% Wb1p3_theory = k2 * dxc * Wc_real / lc * [1/2; ones(Nc - 1, 1); 1/2] + k1 * Wclcp + k3 * Wb_real(:, Nb + 1) / lb + k4 * Wblbp;
% Wb1p3_theory ./ Wblbp3_theory / lb^2
% plot(t_real / lc * sqrt(Hc / mc), Wc_real(:, end) / lc);
% hold on
% plot(t_real / lc * sqrt(Hc / mc), Wc_real(:, 51) / lc);
% hold on
% plot(t_real / lc * sqrt(Hc / mc), Wb_real(:, 51) / lb);
% hold off