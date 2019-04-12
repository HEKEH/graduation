
Configs;

domega = 0.001;
omega = domega : domega: 1;
ret = determinant(omega, beta1, theta, H, lambda, P, k1, k2, k3, k4);
plot(omega, ret)
% 
% [~, min_idx] = min(ret);
% [~, max_idx] = max(ret);
% min_omega = omega(min_idx);
% max_omega = omega(max_idx);

% options=optimset('tolfun',1e-10);
% omega_ans = fsolve(@(omega)determinant(omega, m, beta1, theta, H, lambda, P), 5, options);
% 
% determinant(omega_ans, m, beta1, theta, H, lambda, P)