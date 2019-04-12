%parpool('local', 4)
run('../Configs.m');
run('../MySettings.m');

%params of cable
Jc = eye(Nc - 1);
Lc = zeros(Nc - 1, Nc + 1);
for i = 1: Nc - 1
    Lc(i, i) = 1;
    Lc(i, i + 1) = -2;
    Lc(i, i + 2) = 1;
end
Lc = Lc / dxc^2;

Bc = ones(Nc - 1, Nc + 1) * dxc;
Bc(:, 1) = 1/2 * Bc(:, 1);
Bc(:, end) = 1/2 * Bc(:, end);

%params of beam
Jb = eye(Nb - 3);

Lb = zeros(Nb - 3, Nb + 1);
for i = 1: Nb - 3
    Lb(i, i + 1) = 1;
    Lb(i, i + 2) = -2;
    Lb(i, i + 3) = 1;
end
Lb = Lb / dxb^2;

Gb = zeros(Nb - 3, Nb + 1);
for i = 1: Nb - 3
    Gb(i, i) = 1;
    Gb(i, i + 1) = -4;
    Gb(i, i + 2) = 6;
    Gb(i, i + 3) = -4;
    Gb(i, i + 4) = 1;
end
Gb = Gb / dxb^4;

%solve
Kc1 = zeros(2 * Nc - 2, 2 * Nc);
Kc1(1: Nc - 1, Nc + 2: 2 * Nc) = Jc;
Kc1(Nc: 2 * Nc - 2, 1: Nc + 1) = Lc - lambda^2 * Bc;
Kc1(Nc: 2 * Nc - 2, Nc + 2: 2 * Nc) = - muc * Jc;
Kc2 = zeros(2 * Nc - 2, Nc - 1);
Kc2(Nc: 2 * Nc - 2, :) = eye(Nc - 1);
Kc3 = -lambda^2 * (H * tan(theta) - 1/2) * Kc2;

Kb1 = zeros(2 * Nb - 6, 2 * Nb - 2);
Kb1(1: Nb - 3, Nb + 2: 2 * Nb - 2) = Jb;
Kb1(Nb - 2: 2 * Nb - 6, 1: Nb + 1) = - beta * Gb - P * Lb;
Kb1(Nb - 2: 2 * Nb - 6, Nb + 2: 2 * Nb - 2) = - mub * Jb;
Kb2 = zeros(2 * Nb - 6, Nb - 3);
Kb2(Nb - 2: 2 * Nb - 6, :) = eye(Nb - 3);

% xcs = (1: Nc - 1) * dxc * lc;
% xbs = (2: Nb - 2) * dxb * lb;
% param_c = lc / Hc;
% param_b = mc * lc / (mb * Hc * cos(theta));
% param_t = lc * sqrt(mc / Hc);
% Fc = @(t)fc(xcs, t * param_t) * param_c;
% Fb = @(t)fb(xbs, t * param_t) * param_b;

Fc = @(t)fc((1: Nc - 1) * dxc, t);%这里的fc是无量纲的
Fb = @(t)fb((2: Nb - 2) * dxb, t, m);%这里的fb是无量纲的

param_c = 1/ (lc / Hc);
param_b = 1/ (mc * lc / (mb * Hc * cos(theta)));
param_t = 1/ (lc * sqrt(mc / Hc));

[K_of_fR, row_vector_fR6] = get_params_of_fR(beta, lambda, H, theta, m, dxc, dxb);
R_bar = zeros(1, 2 * Nc + 2 * Nb - 8);
R_bar_p = zeros(1, 2 * Nc + 2 * Nb - 8);
if ~exist('数据','dir')
    mkdir('数据');
end
dRdt = @(t, R_bar, R_bar_p) dRbar_dt(t, R_bar, R_bar_p, Nc, Nb, K_of_fR, row_vector_fR6, H, lambda, Lc, Bc, dxc, theta, Kc1, Kc2, Kc3, Kb1, Kb2, Fc, Fb);
for r = 1: rounds
    t_inteval = (t_per_round * (r - 1)): dt: (t_per_round * r);
    [~, R_bar_p] = decic(dRdt,t_per_round * (r - 1), R_bar(end, :), ones(1, size(R_bar, 2)), zeros(1, size(R_bar, 2)), []); % 求R_bar_p
    [t, R_bar] = ode15i(dRdt, t_inteval, R_bar(end, :), R_bar_p);
    [Rc, Rb] = fR(R_bar', Nc, Nb, K_of_fR, row_vector_fR6);

    Wc = Rc(1: Nc + 1, :);
    Vc = Rc(Nc + 2: 2 * Nc, :);
    Wb = Rb(1: Nb + 1, :);
    Vb = Rb(Nb + 2: 2 * Nb - 2, :);

%     t_real = lc * sqrt(mc / Hc) * t';
%     Wc_real = lc * Wc;
%     Vc_real = Vc / sqrt(mc / Hc);
%     Wb_real = lb * Wb;
%     Vb_real = lb * Vb / (lc * sqrt(mc / Hc));
    writeDatas(t', Wc, Vc, Wb, Vb, r, dNt);
end
