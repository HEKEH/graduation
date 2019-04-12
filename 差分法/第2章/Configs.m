%common
g = 9.81;
theta = 30 / 180 * pi;

%cable: real
Ec = 200 * 1e9;
Ac = 3.14 * 10^(-2);
rhoc = 7850; %kg/m3
mc = rhoc * Ac;

Hc = round(1860 * 10^6 * Ac * 0.35 / 1e7) * 1e7; %ÔÝ¶¨
lc = 100;
%cc = 1 / 100; %ÔÝ¶¨×èÄá
zc = @(xc)1/2 * (mc * g * cos(theta) / Hc .* xc .* (lc - xc));
zc_prime = @(xc)(-mc * g / Hc * cos(theta) .* (xc - lc / 2));
s_prime_3 = @(xc)(sqrt(1 + zc_prime(xc) .^ 2)) .^ 3;
Le = integral(s_prime_3, 0, lc);

%cable: no dimension
muc = 1e-3;
cc = muc / lc * mc / sqrt(mc / Hc);
H = Hc / (mc * g * lc * cos(theta));
lambda = sqrt((mc * g * cos(theta) * lc / Hc)^2 * (Ec * Ac * lc) / (Hc * Le));

%beam: real
Eb = 3.45e10;%N/m2
%Ab = (100 * 50 - 96 * 46) * 1e-6; %m2

% K = 100;
% Ib = K * Ec * Ac^2 / Eb; %m4
Ib = 1;


lb = lc * cos(theta); %m
mb = 4 * 10^4; %kg/m
%cb = 1 / 100; %ÔÝ¶¨×èÄá

%beam: no dimension
mub = 1e-3;
cb = mub / lc * mb / sqrt(mc / Hc);
%mub = cb * lc / mb * sqrt(mc / Hc);

%common: no dimension
m = mc / mb;
P = mc / mb / cos(theta);
beta = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);