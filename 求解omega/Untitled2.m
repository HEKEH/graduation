Ib = K * Ec * Ac^2 / Eb; %m4
beta1 = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);
k1 =  m / (beta1 * cos(theta));
k2 =  m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta));
k3 =  m * (- 2 * H * tan(theta) + 1)^2 * lambda^2 * cos(theta) / (4 * beta1);
k4 =  -m / (beta1 * cos(theta));
func = @(omegas) determinant(omegas, beta1, theta, H, lambda, P, k1, k2, k3, k4);
func2 = @(omegas) determinant2(omegas, beta1, theta, H, lambda, P, k1, k2, k3, k4);

domega = 1e-4;
omegas = domega: domega: uplimit;
eqn_ans = func(omegas);
eqn_ans2 = func2(omegas);
eqn_ans ./ eqn_ans2;

% low_idxs = (eqn_ans(1: end - 1) .* eqn_ans(2: end)) <= 0;
% lows = omegas(low_idxs);
% highs = lows + domega;
% 
% for i = length(lows): -1: 1
%     if func(highs(i)) > func(lows(i))
%         if (func(lows(i)) < func(lows(i) - 1e-15)) && (func(highs(i)) > func(highs(i) + 1e-15))
%             lows(i) = [];
%             highs(i) = [];
%         end
%     else
%         if (func(lows(i)) > func(lows(i) - 1e-15)) && (func(highs(i)) < func(highs(i) + 1e-15))
%             lows(i) = [];
%             highs(i) = [];
%         end
%     end
% end
% 
% omegas_num = length(lows);
% omega_ans = zeros(1, omegas_num);
% for i = 1: omegas_num
%     low = lows(i);
%     high = highs(i);
%     incline = func(high) >= func(low);
%     mid = (low + high) / 2;
%     limit = 1e-5;
%     rounds = 0;
%     while (abs(func(mid)) > limit)
%         rounds = rounds + 1;
%         if rounds > 200
%             low = lows(i);
%             high = highs(i);
%             mid = (low + high) / 2;
%             rounds = 0;
%             limit = limit * 10;
%         end
%         if func(mid) > 0
%             if incline
%                 high = mid;
%             else
%                 low = mid;
%             end
%         else
%             if incline
%                 low = mid;
%             else
%                 high = mid;
%             end
%         end
%         mid = (low + high) / 2;
%     end
%     omega_ans(i) = mid;
%     limit
% end
    %determinant(omega_ans, beta1, theta, H, lambda, P, k1, k2, k3, k4)