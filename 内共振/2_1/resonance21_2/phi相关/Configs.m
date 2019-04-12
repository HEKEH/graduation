%digits(200);
g = 9.81; %m/s2
theta = 30 / 180 * pi;

%cable: real
Ec = 200 * 1e9;
lc = 100;
Ac = 3.14 * 10^(-2);
mc = 7850 * Ac;
Hc = round(1860 * 10^6 * Ac * 0.35 / 1e7) * 1e7; %ÔÝ¶¨
%Hc = 1;

%beam: real
Eb = 3.45e10;%N/m2
%Ib = K * Ec * Ac^2 / Eb; %m4
lb = lc * cos(theta); %m
mb = 4 * 10^4; %kg/m
zc = @(xc)1/2 * (mc * g * cos(theta) / Hc .* xc .* (lc - xc));
zc_prime = @(xc)(-mc * g / Hc * cos(theta) .* (xc - lc / 2));
s_prime_3 = @(xc)(sqrt(1 + zc_prime(xc) .^ 2)) .^ 3;
Le = integral(s_prime_3, 0, lc);

%cable: no dimension
H = Hc / (mc * g * lc * cos(theta));
lambda = sqrt((mc * g * cos(theta) * lc / Hc)^2 * (Ec * Ac * lc) / (Hc * Le));

%beam: no dimension

%common: no dimension
m = mc / mb;
P = mc / mb / cos(theta);
%beta1 = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);