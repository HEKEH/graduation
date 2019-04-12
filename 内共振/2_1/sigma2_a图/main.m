clear;
choice = 'n';
run(['./',choice, '/get_p.m']);
sigma2_range = -2e-2: 1e-5: 2e-2;
sigma2_len = length(sigma2_range);
a_cell = cell(2, sigma2_len);%第一行m，第二行n
gamma1_angle = cell(1, sigma2_len);

if choice == 'm'
    get_a_n = @(a_m, sigma2) [-((a_m.^2.*pn(2)+4.*sigma1.*omega_n+4.*sigma2.*omega_n-2.*sqrt(a_m.^2.*(m.*J_n+pn(4).*sigma1).^2-4.*((-1+pn(1)).*mu_b-pn(1).*mu_c).^2.*omega_n.^2))./pn(3));...
              -((a_m.^2.*pn(2)+4.*sigma1.*omega_n+4.*sigma2.*omega_n+2.*sqrt(a_m.^2.*(m.*J_n+pn(4).*sigma1).^2-4.*((-1+pn(1)).*mu_b-pn(1).*mu_c).^2.*omega_n.^2))./pn(3))];
    get_ret = @(a_m, a_n, sigma2) -a_m.^2.*(pm(5)+pm(6)).^2.*(m.*J_n+pn(4).*sigma1).^2+1./4.*(a_m.^2.*(((-(((-1)+pm(1)))).*mu_b+pm(1).*mu_c)).*((m.*J_n+pn(4).*sigma1)).*omega_m+a_n.^2.*(((-(((-1)+pn(1)))).*mu_b+pn(1).*mu_c)).*((m.*J_m+pm(4).*sigma1)).*omega_n).^2+1./64.*(a_n.^2.*((a_m.^2.*pn(2)+a_n.^2.*pn(3))).*((m.*J_m+pm(4).*sigma1))-a_m.^2.*((a_m.^2.*pm(2)+a_n.^2.*pm(3))).*((m.*J_n+pn(4).*sigma1))-8.*a_m.^2.*(m.*J_n+pn(4).*sigma1).*sigma2.*omega_m+4.*a_n.^2.*(m.*J_m+pm(4).*sigma1).*(sigma1+sigma2).*omega_n).^2;
    get_cos_gamma1 = @(a_m, a_n, sigma2) -((1./2.*(a_m.^2.*pn(2)+a_n.^2.*pn(3))+2.*omega_n(sigma1+sigma2))./a_m(m.*J_n+pn(4).*sigma1));
    get_sin_gamma1 = @(a_m, a_n, sigma2) (2*((1-pn(1)).*mu_b+pn(1).*mu_c).*omega_n)./(a_m.*(m.*J_n+pn(4).*sigma1));
else
    get_a_n = @(a_m, sigma2)sqrt([(a_m.^2.*pm(3).*((a_m.^2.*pm(2)-8.*((sigma1-2.*sigma2)).*omega_m))+2.*sqrt(a_m.^6.*pm(2).^2.*(m.*J_m+pm(4).*sigma1).^2+16.*a_m.^2.*(m.*J_m+pm(4).*sigma1).^2.*(((-1+pm(1)).*mu_b-pm(1).*mu_c).^2+4.*(sigma1-2.*sigma2).^2).*omega_m.^2+4.*a_m.^4.*omega_m.*(-4.*pm(2).*(m.*J_m+pm(4).*sigma1).^2.*(sigma1-2.*sigma2)-pm(3).^2.*((((-1)+pm(1))).*mu_b-pm(1).*mu_c).^2.*omega_m)))./((2.*m.*J_m-a_m.*pm(3)+2.*pm(4).*sigma1).*(2.*m.*J_m+a_m.*pm(3)+2.*pm(4).*sigma1));...
              (a_m.^2.*pm(3).*((a_m.^2.*pm(2)-8.*((sigma1-2.*sigma2)).*omega_m))-2.*sqrt(a_m.^6.*pm(2).^2.*(m.*J_m+pm(4).*sigma1).^2+16.*a_m.^2.*(m.*J_m+pm(4).*sigma1).^2.*(((-1+pm(1)).*mu_b-pm(1).*mu_c).^2+4.*(sigma1-2.*sigma2).^2).*omega_m.^2+4.*a_m.^4.*omega_m.*(-4.*pm(2).*(m.*J_m+pm(4).*sigma1).^2.*(sigma1-2.*sigma2)-pm(3).^2.*((((-1)+pm(1))).*mu_b-pm(1).*mu_c).^2.*omega_m)))./((2.*m.*J_m-a_m.*pm(3)+2.*pm(4).*sigma1).*(2.*m.*J_m+a_m.*pm(3)+2.*pm(4).*sigma1))]);
    get_ret = @(a_m, a_n, sigma2) -a_n.^2.*(pn(5)+pn(6)).^2.*(m.*J_m+pm(4).*sigma1).^2+1./4.*(a_m.^2.*(((((-1)+pm(1))).*mu_b-pm(1).*mu_c)).*((m.*J_n+pn(4).*sigma1)).*omega_m-a_n.^2.*(((((-1)+pn(1))).*mu_b-pn(1).*mu_c)).*((m.*J_m+pm(4).*sigma1)).*omega_n).^2+1./64.*(a_m.^2.*((m.*J_n+pn(4).*sigma1)).*((a_m.^2.*pm(2)+a_n.^2.*pm(3)-8.*((sigma1-2.*sigma2)).*omega_m))-a_n.^2.*((m.*J_m+pm(4).*sigma1)).*((a_m.^2.*pn(2)+a_n.^2.*pn(3)+8.*sigma2.*omega_n))).^2;
    get_cos_gamma1 = @(a_m, a_n, sigma2)-((a_m.*(1./2.*(a_m.^2.*pm(2)+a_n.^2.*pm(3))+4.*(-sigma1+2.*sigma2).*omega_m))./(a_n.^2.*((m.*J_m+pm(4).*sigma1))));
    get_sin_gamma1 = @(a_m, a_n, sigma2)-((2*((1-pm(1)).*mu_b+pm(1).*mu_c).*a_m.*omega_m)./(a_n.^2.*((m.*J_m+pm(4).*sigma1))));
end

valid_idx = [];
for sidx = 1: sigma2_len
    sigma2 = sigma2_range(sidx);
    [a_m,a_n] = get_a(get_a_n, get_ret, sigma2, @get_ret_2);
    if ~isempty(a_m)
        valid_idx = [valid_idx, sidx];
        a_cell{1, sidx} = a_m;
        a_cell{2, sidx} = a_n;
        
        tmp = get_cos_gamma1(a_m, a_n, sigma2);
        tmp(tmp < -1) = -1;
        tmp(tmp > 1) = 1;
        angle = acos(tmp);
        nn_idx = get_sin_gamma1(a_m, a_n, sigma2) < 0;
        angle(nn_idx) =  - angle(nn_idx);
        gamma1_angle{1, sidx} = angle;
    end
end
if isempty(valid_idx)
    disp('无有效数据');
else
    sigma2_range = sigma2_range(valid_idx);
    a_cell = a_cell(:, valid_idx);
    sigma2_len = length(sigma2_range);

    left = zeros(2, 1);
    %PL = zeros(1, 2);
    %lgd = cell(1, 2);
    lgd = {'a_m', 'a_n'};
    max_idx = zeros(2, 1);
    colors = {'r', [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250],...
        [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
    list = cell(2,1);
    
    %画点
    hold off
    for idx = 1: 2
        figure(idx)
        for sidx = 1: sigma2_len
            plot(sigma2_range(sidx), a_cell{idx, sidx}, '.', 'color', colors{idx})
            hold on
        end
    end
    hold off
    %画点
    

    
%     for idx = 1: 2
%         [list{idx}, max_size, start_idx, end_idx, left(idx)] = cell_to_list(sigma2_len, a_cell, idx);
%         if isempty(list{idx})
%             break;
%         end
%         figure(idx)
%         linewidth = 2;
%         if max_size == 1
%             [~, max_idx(idx)] = max(list{idx});
%             PL = plot(sigma2_range, list{idx}, 'color', colors{idx} ,'LineWidth', linewidth);
%         elseif max_size == 2
%             [~, max_idx(idx)] = max(max(list{idx}));
%             PL = plot(sigma2_range, list{idx}, '.', 'color', colors{idx} ,'LineWidth', linewidth);
%         elseif max_size == 3
%             [~, max_idx(idx)] = max(max(list{idx}));
%             PL = plot(sigma2_range(1: end_idx), list{idx}(1, 1: end_idx), 'color', colors{idx} ,'LineWidth', linewidth);
%             hold on
%             plot(sigma2_range(start_idx: end_idx), list{idx}(2, start_idx: end_idx), ':', 'color', colors{idx}, 'LineWidth', linewidth);
%             hold on
%             plot(sigma2_range(start_idx: sigma2_len), list{idx}(3, start_idx: sigma2_len), 'color', colors{idx} ,'LineWidth', linewidth);
%         end
%         hold off
%         legend(PL, lgd{idx})
%         set(gca,'FontSize',14)
%         set (gcf,'Position',[50,50,1000,750], 'color','w') 
%         set(gca,'Position',[0.075 0.1 0.85 0.85]);
%         xlabel('\it\sigma_2', 'FontSize', 15)
%         ylabel(['\it', lgd{idx}], 'FontSize', 15)
%         ylim([0, Inf]);
%     end
%     %暂时,没考虑多重解
%     a_m = list{1};
%     a_n = list{2};
%     
%     [gamma1_list, max_size, start_idx, end_idx, ~] = cell_to_list(sigma2_len, gamma1_angle,1); %idx, left(idx));
%     
%     a_bar(1,:) = sqrt(a_m.^2+(cos(gamma1_list).*a_m.*a_n.^2.*((cos(theta).^3.*integral(@(xb)phi_b(xb,m_idx).*psi_bks{1}(xb), 0, 1)+m.*integral(@(xc)phi_c(xc,m_idx).*psi_cks{1}(xc), 0, 1))))+(1./4.*a_n.^4.*((cos(theta).^3.*integral(@(xb)psi_bks{1}(xb).^2, 0, 1)+m.*integral(@(xc)psi_cks{1}(xc).^2, 0, 1)))));
%     a_bar(2,:) = sqrt(a_n.^2+(cos(gamma1_list).*a_m.*a_n.^2.*((cos(theta).^3.*integral(@(xb)phi_b(xb,n_idx).*psi_bks{2}(xb), 0, 1)+m.*integral(@(xc)phi_c(xc,n_idx).*psi_cks{2}(xc), 0, 1))))+(1./4.*a_m.^2.*a_n.^2.*((cos(theta).^3.*integral(@(xb)psi_bks{2}(xb).^2, 0, 1)+m.*integral(@(xc)psi_cks{2}(xc).^2, 0, 1)))));
%     
%     for idx = 1: 2
%         figure(idx + 2)
%         linewidth = 2;
%         PL = plot(sigma2_range, a_bar(idx, :), 'color', colors{idx} ,'LineWidth', linewidth);
%         hold off
%         legend(PL, lgd{idx})
%         set(gca,'FontSize',14)
%         set (gcf,'Position',[50,50,1000,750], 'color','w') 
%         set(gca,'Position',[0.075 0.1 0.85 0.85]);
%         xlabel('\it\sigma_2', 'FontSize', 15)
%         ylabel(['\it', lgd{idx}], 'FontSize', 15)
%         ylim([0, Inf]);
%     end
end
