if ~exist('a_cell', 'var')
    if ~exist('data.mat', 'file')
        Main;
        save data.mat
    else 
        load('data.mat')
    end
end
idx = 1;

sigma_idx = max_idxs(idx);
this_sigma = sigma_range(sigma_idx);
[this_a, temp] = max(a_cell{idx, sigma_idx});
this_gamma = gamma_angle{idx, sigma_idx}(temp);


this_OMEGA = omegas(idx) + this_sigma;
wc1 = @(xc, t) this_a * 2 * cos(this_OMEGA * t - this_gamma) .* phi_c(xc, idx);
wc2 = @(xc, t) this_a^2 * 2 * (cos(2 * (this_OMEGA * t - this_gamma)) .* psi_c2(xc, idx) + psi_c3(xc, idx));

wb1 = @(xb, t) this_a * 2 * cos(this_OMEGA * t - this_gamma) .* phi_b(xb, idx);
wb2 = @(xb, t) this_a^2 * 2 * (cos(2 * (this_OMEGA * t - this_gamma)) .* psi_b2(xb, idx) + psi_b3(xb, idx));

wc = @(xc, t) wc1(xc, t) + wc2(xc, t);
wb = @(xb, t) wb1(xb, t) + wb2(xb, t);

enlarge_times = 2e1;%放大倍数
this_x = 0: 1e-2: 1;
for t = 0: 0.2: 1000
    hold off
    plot(lb * this_x + enlarge_times * (lb * wb(1, t) * this_x * sin(theta) * cos(theta) - lc * wc(this_x, t) * sin(theta)), ...
        lb * (1 - this_x) * tan(theta) - enlarge_times * (lc * wc(this_x, t) * cos(theta) + lb * wb(1, t) * this_x * sin(theta) * sin(theta)),...
        'linewidth', 1.5);
    hold on
    plot(lb * this_x, - enlarge_times * lb * wb(this_x, t), 'linewidth', 3);
    axis equal
    axis([0, ceil(lb / 10) * 10, -10, Inf])
    title(['实际时间: ', num2str(t * lc / sqrt(Hc / mc)), 's 无量纲时间: ', num2str(t), ' (振动放大', num2str(enlarge_times), '倍)']);
    drawnow;
end
plot(0: 0.2: 1000,wc(0.5, 0: 0.2: 1000))