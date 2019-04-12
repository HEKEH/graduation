idx = 200;
K = Ks(idx);
Ib = K * Ec * Ac^2 / Eb; %m4
beta1 = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);
k1 =  m / (beta1 * cos(theta));
k2 =  m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta));
k3 =  m * (- 2 * H * tan(theta) + 1)^2 * lambda^2 * cos(theta) / (4 * beta1);
k4 =  -m / (beta1 * cos(theta));
func2 = @(omegas) determinant2(omegas, beta1, theta, H, lambda, P, k1, k2, k3, k4);

domega = 1e-5;
omegas = domega: domega: 18.81;
eqn_ans = func2(omegas);
plot(omegas, eqn_ans);

% mids = ret{idx};
% domega = 1e-6;
% for i = 1: 3
%     figure(i)
%     omegas = max(domega,mids(i) - 1e2 * domega): domega: mids(i) + 1e2 * domega;
%     eqn_ans = func(omegas);
%     plot(omegas, eqn_ans,'.');
%     hold on 
%     plot(omegas, zeros(1, length(omegas)));
%     hold on 
%     plot(mids(i), func(mids(i)),'*')
%     hold off
% end


