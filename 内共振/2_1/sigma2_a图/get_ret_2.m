function ret = get_ret_2(a_m, idx, sigma2, get_ret, get_a_n)
    a_n_all = get_a_n(a_m, sigma2);
    a_n = a_n_all(idx);
    ret = get_ret(a_m, a_n, sigma2);
end