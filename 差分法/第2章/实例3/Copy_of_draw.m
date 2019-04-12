%parpool('local', 4)
Configs;
MySettings;
folder = '数据';
for r = 1: rounds
    t_real = csvread([folder, '\t',num2str(r),'.csv']); 
    Wc_real = csvread([folder, '\Wc',num2str(r),'.csv']); 
    Vc_real = csvread([folder, '\Vc',num2str(r),'.csv']); 
    Wb_real = csvread([folder, '\Wb',num2str(r),'.csv']); 
    Vb_real = csvread([folder, '\Vb',num2str(r),'.csv']); 
    enlarge_times = 5e2;%放大倍数
    for i = 1: 30: length(t_real)
        hold off
        plot(lb * (0: dxc: 1) + enlarge_times * (Wb_real(i, Nb) * (0: dxc: 1) * sin(theta) * cos(theta) - Wc_real(i, :) * sin(theta)), ...
            lb * (1: -dxc: 0) * tan(theta) - enlarge_times * (Wc_real(i, :) * cos(theta) + Wb_real(i, Nb) * (0: dxc: 1) * sin(theta) * sin(theta)),...
            'linewidth', 1.5);
        hold on
        plot(lb * (0: dxb: 1), - enlarge_times * Wb_real(i, :), 'linewidth', 3);
        axis([-lb * 0.1, lb * 1.1, -0.2 * lb, lb])
        title(['t = ', num2str(t_real(i)), 's (振动放大', num2str(enlarge_times), '倍)']);
        drawnow;
    end
end