Configs;
pow = 0: 0.0005: 4;
Ks = 10 .^ pow;
len = length(Ks);
uplimit = 20;
ret = cell(1, len);
for i = 1: len
    K = Ks(i);
    temp = get_omegas(K, lb, Eb, lc, Ec, Ac, Hc, mb, mc, m, theta, H, lambda, P, uplimit);
    ret{i} = temp(1: end);
end
% for i = 1: size(ret, 1)
%     for j = 1: size(ret, 2)
%         min_d = 100;
%         if i ~= 1
%             min_d = min(min_d, abs(ret(i, j) - ret(i, j - 1));
%         end
%     end
% end

for i = 1: len
    plot(Ks(i), ret{i} / pi, 'r.','MarkerSize',1);
    hold on
end
%legend(P, lgd, 'location','northwest', 'FontSize', 12);
set(gca,'xscale','log');
set(gca,'FontSize',14)
set (gcf,'Position',[50,50,1000,750], 'color','w') 
set(gca,'Position',[0.075 0.1 0.85 0.85]);
xlabel('¸Õ¶È±ÈK', 'FontSize', 15)
ylabel('ÆµÂÊ\omega/\pi', 'FontSize', 15)
axis([0, max(Ks), 0, uplimit / pi]);

hold off

two_to_one_omegas = [];
for idx = 1:length(ret)
    temp = ret{idx};
    if abs(temp(2) / temp(1) - 2) < 1e-3
        two_to_one_omegas = [two_to_one_omegas, [Ks(idx); temp(1); temp(2)]];
    end
end
