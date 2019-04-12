%parpool('local', 4)
clear;
if ~exist('data.mat', 'file')
    run('../Configs');
    run('../MySettings');
    
    folder = '数据';
    interval = 5;

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
    save data.mat
else
    if ~exist('Ac', 'var')
        load('data.mat')
    end
end

enlarge_times = 5;%放大倍数
Wc_real = Wc * lc;
Wb_real = Wb * lb;
for i = 1: 1: length(t)
    hold off
    plot(lb * (0: dxc: 1) + enlarge_times * (Wb_real(i, Nb) * (0: dxc: 1) * sin(theta) * cos(theta) - Wc_real(i, :) * sin(theta)), ...
        lb * (1: -dxc: 0) * tan(theta) - enlarge_times * (Wc_real(i, :) * cos(theta) + Wb_real(i, Nb) * (0: dxc: 1) * sin(theta) * sin(theta)),...
        'linewidth', 1.5);
    hold on
    plot(lb * (0: dxb: 1), - enlarge_times * Wb_real(i, :), 'linewidth', 3);
    axis equal
    axis([0, ceil(lb / 10) * 10, -10, Inf])
    title(['无量纲时间: ', num2str(t(i)), ' 实际时间: ', num2str(t(i) * lc / sqrt(Hc / mc)), 's (振动放大', num2str(enlarge_times), '倍)']);
    %xlim([0, ceil(lb / 10) * 10]);
    drawnow;
end

% hc = Ec * Ac / Le * ((sin(theta) - mc * g * lc * cos(theta)^2 / 2 / Hc) * Wb(:, Nb + 1) + mc * g * cos(theta) / Hc * Wc * [1/2; ones(Nc - 1, 1); 1/2] * dxc * lc);
% tan_phi = mc * g * lc * cos(theta) / 2 / Hc;
% Wblbp = (3 * Wb(:, Nb + 1) - 4 * Wb(:, Nb) + Wb(:, Nb - 1)) / (2 * dxb * lb);
% Wclcp = (3 * Wc(:, Nc + 1) - 4 * Wc(:, Nc) + Wc(:, Nc - 1)) / (2 * dxc * lc);
% Wblbp3 = (3 * Wb(:, Nb - 3) - 14 * Wb(:, Nb - 2) + 24 * Wb(:, Nb - 1) - 18 * Wb(:, Nb) + 5 * Wb(:, Nb + 1)) / (2 * (dxb * lb)^3);
% Wblbp3_theory = (hc * (sin(theta) - cos(theta) * tan_phi) - Hc * Wblbp * cos(theta) + Hc * Wclcp * cos(theta)) / Eb / Ib;
% Wblbp3 ./ Wblbp3_theory
% Wb1p3_theory = k2 * dxc * Wc / lc * [1/2; ones(Nc - 1, 1); 1/2] + k1 * Wclcp + k3 * Wb(:, Nb + 1) / lb + k4 * Wblbp;
% Wb1p3_theory ./ Wblbp3_theory / lb^2
plot(t, Wc(:, end));
hold on
plot(t, Wc(:, 51));
hold on
plot(t, Wb(:, 51));
hold off