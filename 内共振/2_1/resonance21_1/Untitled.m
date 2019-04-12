% syms omega_n x
% a = sqrt(omega_n);
% arr = [cos(a * x) sin(a * x) cosh(a * x) sinh(a * x)];
% mtx = [subs(arr, x, 0); subs(arr, x, 1);  subs(diff(arr, x, 2), x, 0); subs(diff(arr, x, 2), x, 1)];
% ret = det(mtx);
% null(double(subs(mtx, omega_n, (2 * pi)^2)))

temp_m = double(eval(matr_A))
temp_b = double(eval(vec_b))
