K = 100;
omegas_num = 8;
omegas_idx = 1: omegas_num;
sigma_range = -2e-2: 5e-6: 2e-2;
sigma_len = length(sigma_range);
all_params;

mu_c = 1e-3;
mu_b = 1e-3;
F_amplitude = 1e-4;
Fc = @(xc) F_amplitude .* ones(size(xc, 1), size(xc, 2));
Fb = @(xb) F_amplitude .* m * ones(size(xb, 1), size(xb, 2));%m

p = get_p(omegas_num, omegas_idx, m, theta, PHI_c3, phi_c, phi_b, Fc, Fb);


a_cell = cell(omegas_num, sigma_len);
gamma_angle = cell(omegas_num, sigma_len);
max_idxs = zeros(omegas_num, 1);
for idx = omegas_idx
    temp1 = 1/ 64 * p(idx, 2) ^ 2;
    temp4 = (p(idx, 3) + p(idx, 4))^2;
    for sidx = 1: sigma_len
        sigma1 = sigma_range(sidx);
        a = get_a(p, temp1, temp4, sigma1, omegas, mu_c, mu_b, idx);
        a_cell{idx, sidx} = a;
        %gamma_angle(idx, sidx) = acot((-p (idx, 2) * a^3 - 2 * omegas(idx) * sigma1 * a) / (omegas(idx) * (p(idx, 1) * mu_c + (1 - p(idx, 1)) * mu_b) * a));
        tmp = (-1/8 * p (idx, 2) * a.^3 - omegas(idx) * sigma1 * a) / (p(idx, 3) + p(idx, 4));
        tmp(tmp < -1) = -1;
        tmp(tmp > 1) = 1;
        angle = acos(tmp);
        n_idx = (1/ 2 * omegas(idx) * (p(idx, 1) * mu_c + (1 - p(idx, 1)) * mu_b) * a) / (p(idx, 3) + p(idx, 4)) < 0;
        angle(n_idx) = -angle(n_idx);
        gamma_angle{idx, sidx} = angle;
    end
end
% 
left = zeros(omegas_num, 1);
PL = zeros(1, omegas_num);
lgd = cell(1, omegas_num);
colors = {'r', [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
figure(1)
linewidth = 2;

p_idxs = [3,4,7];
for idx = p_idxs
    [list, max_size, start_idx, end_idx, left(idx)] = cell_to_list(sigma_len, a_cell, idx); %
    if max_size == 1
        [~, max_idxs(idx)] = max(list);
        PL(idx) = plot(sigma_range, list, 'color', colors{idx} ,'LineWidth', linewidth);
    else
        [~, max_idxs(idx)] = max(max(list));
        PL(idx) = plot(sigma_range(1: end_idx), list(1, 1: end_idx), 'color', colors{idx} ,'LineWidth', linewidth);
        hold on
        plot(sigma_range(start_idx: end_idx), list(2, start_idx: end_idx), ':', 'color', colors{idx}, 'LineWidth', linewidth);
        hold on
        plot(sigma_range(start_idx: sigma_len), list(3, start_idx: sigma_len), 'color', colors{idx} ,'LineWidth', linewidth);
    end
    hold on
    plot(0, list(size(list, 1), ceil(length(list) / 2)), 'r.', 'MarkerSize', 15)
    hold on
    lgd{idx} = ['µÚ', num2str(idx), '½×'];
end
hold off
legend(PL(p_idxs), lgd{p_idxs})
set(gca,'FontSize',14)
set (gcf,'Position',[50,50,1000,750], 'color','w') 
set(gca,'Position',[0.075 0.1 0.85 0.85]);
xlabel('\fontname{roman}\sigma\epsilon^2', 'FontSize', 15)
ylabel('\fontname{roman}\it a¡Á\epsilon', 'FontSize', 15)

% figure(2)
% for idx = omegas_idx
%     [list, max_size, start_idx, end_idx, ~] = cell_to_list(sigma_len, gamma_angle, idx, left(idx));
%     if max_size == 1
%         PL(idx) = plot(sigma_range, list, 'color', colors{idx}, 'LineWidth', linewidth);
%     else
%         PL(idx) = plot(sigma_range(1: end_idx), list(1, 1: end_idx), 'color', colors{idx}, 'LineWidth', linewidth);
%         hold on
%         plot(sigma_range(start_idx: end_idx), list(2, start_idx: end_idx), ':', 'color', colors{idx}, 'LineWidth', linewidth);
%         hold on
%         plot(sigma_range(start_idx: sigma_len), list(3, start_idx: sigma_len), 'color', colors{idx}, 'LineWidth', linewidth);
%     end
%     hold on
% end
% hold off
% legend(PL, lgd, 'location', 'northwest', 'FontSize', 12);
% set(gca,'FontSize',14)
% set(gca,'YTick',-pi: pi/2: pi)
% set(gca,'ytickLabel',{'-¦Ð','-¦Ð/2','0','¦Ð/2','¦Ð'})%'-¦Ð','-¦Ð/2',
% ylim([-pi pi])
% set (gcf,'Position',[50,50,1000,750], 'color','w') 
% set(gca,'Position',[0.075 0.1 0.85 0.85]);
% xlabel('\sigma\epsilon^2', 'FontSize', 15)
% ylabel('\it \gamma', 'FontSize', 15)
% 
% sigma_max = sigma_range(max_idxs);