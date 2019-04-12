function [a_m_ans, a_n_ans] = get_a(get_a_n, get_ret, sigma2, get_ret_2)
    d_a_m = 1e-8;
    max_a_m = 1e-3;
    a_m_all = d_a_m: d_a_m: max_a_m;
    a_n_all = get_a_n(a_m_all, sigma2);
    reasonable_idx = (imag(a_n_all) == 0);
    a_m_ans = [];
    a_n_ans = [];
    for idx = 1: 2
        a_m = a_m_all(reasonable_idx(idx,:));
        if ~isempty(a_m)
            get_ret_temp = @(a_m)get_ret_2(a_m, idx, sigma2, get_ret, get_a_n);
            a_n = a_n_all(idx, reasonable_idx(idx,:));
            ret = get_ret(a_m, a_n, sigma2);
            low_idxs = (ret(1: end - 1) .* ret(2: end)) <= 0;
            lows = a_m(low_idxs);
            highs = a_m(find(low_idxs) + 1);
            temp_idx = (abs(highs - lows - d_a_m) > 1/2 * abs(d_a_m));
            highs(temp_idx) = [];
            lows(temp_idx) = [];
            a_len = length(lows);
            for i = 1: a_len
                low = lows(i);
                high = highs(i);
                incline = get_ret_temp(high) > get_ret_temp(low);
                mid = (low + high) / 2;
                limit = 1e-20;
                rounds = 0;
                while (abs(get_ret_temp(mid)) > limit)
                    rounds = rounds + 1;
                    if rounds > 50
                        low = lows(i);
                        high = highs(i);
                        mid = (low + high) / 2;
                        rounds = 0;
                        limit = limit * 10;
                    end
                    if get_ret_temp(mid) > 0
                        if incline
                            high = mid;
                        else
                            low = mid;
                        end
                    else
                        if incline
                            low = mid;
                        else
                            high = mid;
                        end
                    end
                    mid = (low + high) / 2;
                end
                a_m_ans = [a_m_ans; mid];
                temp = get_a_n(mid, sigma2);
                a_n_ans = [a_n_ans; temp(idx)];
            end
        end
    end
end
