K = 100;
omegas_num = 2;
omegas_idx = 2;
sigma_range = -2e-2: 5e-6: 2e-2;
sigma_len = length(sigma_range);
all_params;

mu_c = 1e-3;
mu_b = 1e-3;
%F_amplitude = 1e-2;
% Fc = @(xc) 1e-3 * phi_c(xc, 1) .* ones(size(xc, 1), size(xc, 2));
% Fb = @(xb) 1e-3 * phi_b(xb, 1) .* ones(size(xb, 1), size(xb, 2));%m

fidxs = 1: 5;
Fc = cell(1, length(fidxs));
Fb = cell(1, length(fidxs));
PL = zeros(1, length(fidxs));
lgd = cell(1, length(fidxs));

F_temp = 5e-6;
Fc{1} = @(xc) 2 * F_temp * phi_c(xc, omegas_idx) .* ones(size(xc, 1), size(xc, 2));
Fb{1} = @(xb) 2 * F_temp  * phi_b(xb, omegas_idx) .* ones(size(xb, 1), size(xb, 2));%m
Fc{2} = @(xc) 4 * F_temp  * phi_c(xc, omegas_idx) .* ones(size(xc, 1), size(xc, 2));
Fb{2} = @(xb) 4 * F_temp  * phi_b(xb, omegas_idx) .* ones(size(xb, 1), size(xb, 2));%m
Fc{3} = @(xc) 6 * F_temp  * phi_c(xc, omegas_idx) .* ones(size(xc, 1), size(xc, 2));
Fb{3} = @(xb) 6 * F_temp  * phi_b(xb, omegas_idx) .* ones(size(xb, 1), size(xb, 2));%m
Fc{4} = @(xc) 8 * F_temp  * phi_c(xc, omegas_idx) .* ones(size(xc, 1), size(xc, 2));
Fb{4} = @(xb) 8 * F_temp  * phi_b(xb, omegas_idx) .* ones(size(xb, 1), size(xb, 2));%m
Fc{5} = @(xc) 10 * F_temp  * phi_c(xc, omegas_idx) .* ones(size(xc, 1), size(xc, 2));
Fb{5} = @(xb) 10 * F_temp  * phi_b(xb, omegas_idx) .* ones(size(xb, 1), size(xb, 2));%m
% Fc{5} = @(xc) 1e-3 * phi_c(xc, omegas_idx) .* ones(size(xc, 1), size(xc, 2));
% Fb{5} = @(xb) 1e-3 * phi_b(xb, omegas_idx) .* ones(size(xb, 1), size(xb, 2));%m


colors = {'r', [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
for fidx = fidxs
    p = get_p(omegas_num, omegas_idx, m, theta, PHI_c3, phi_c, phi_b, Fc{fidx}, Fb{fidx});

    a_cell = cell(omegas_num, sigma_len);
    gamma_angle = cell(omegas_num, sigma_len);
    max_idxs = zeros(omegas_num, 1);
    idx = omegas_idx; %%%
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

    % 
    left = zeros(omegas_num, 1);

    %figure(1)
    linewidth = 2;
    [list, max_size, start_idx, end_idx, left(idx)] = cell_to_list(sigma_len, a_cell, idx); %
    if max_size == 1
        [~, max_idxs(idx)] = max(list);
        PL(fidx) = plot(sigma_range, list, 'color', colors{fidx} ,'LineWidth', linewidth);
    else
        [~, max_idxs(idx)] = max(max(list));
        PL(fidx) = plot(sigma_range(1: end_idx), list(1, 1: end_idx), 'color', colors{fidx} ,'LineWidth', linewidth);
        hold on
        plot(sigma_range(start_idx: end_idx), list(2, start_idx: end_idx), ':', 'color', colors{fidx}, 'LineWidth', linewidth);
        hold on
        plot(sigma_range(start_idx: sigma_len), list(3, start_idx: sigma_len), 'color', colors{fidx} ,'LineWidth', linewidth);
    end
    hold on
    lgd{fidx} = ['¹¤¿ö', num2str(fidx)];
end
legend(PL, lgd)
set(gca,'FontSize',14)
set (gcf,'Position',[50,50,1000,750], 'color','w') 
set(gca,'Position',[0.075 0.1 0.85 0.85]);
xlabel('\fontname{roman}\sigma\epsilon^2', 'FontSize', 15)
ylabel('\fontname{roman}\it a¡Á\epsilon', 'FontSize', 15)
hold off

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