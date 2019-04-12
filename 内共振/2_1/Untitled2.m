syms Ups_c11 xc omega_m Ups_c12 Ups_c13 Ups_c14 omega_n Ups_c15 xb Ups_b11 Ups_b12 Ups_b13 Ups_b14 alpha_b1 beta_b1 k1 k2 k3 k4 Gamma1 Gamma2 C_idx
% f = Ups_c11 * xc * cos(omega_m * xc) + Ups_c12 * xc * sin(omega_m * xc)...
%     + Ups_c13 * cos(omega_n * xc) + Ups_c14 * sin(omega_n * xc) + Ups_c15;
% diff(f, xc)
% diff(f, xc, 2)
% int(f,'xc',0, 1)
% g = Ups_b11 * xb .* cos(alpha_b1 * xb) + Ups_b12 * xb .* sin(alpha_b1 * xb)...
%     + Ups_b13 * xb .* cosh(beta_b1 * xb) + Ups_b14 * xb .* sinh(beta_b1 * xb);
% diff(g, xb)
% diff(g, xb, 2)
% diff(g, xb, 3)

c = [cos(omega_m * xc)+Gamma1, sin(omega_m * xc)+Gamma2, 0, 0, 0, 0];
b = [0, 0, cos(alpha_b1 * xb), sin(alpha_b1 * xb), cosh(beta_b1 * xb), sinh(beta_b1 * xb)];
cp = diff(c, xc);
cint = int(c,'xc',0, 1);
bp3 = diff(b, xb, 3);
bp = diff(b, xb, 1);

ret = k1 * cp + k2 * cint + k3 * b + k4 * bp - bp3;
subs(ret, {xb, xc}, {1, 1})