function omega_ans = get_omegas(func, omega_nums)
    domega = 5e-4;
    omegas = domega: domega: 100;
    eqn_ans = func(omegas);
    low_idxs = (eqn_ans(1: end - 1) .* eqn_ans(2: end)) <= 0;
    lows = omegas(low_idxs);
    highs = lows + domega;
    for i = length(lows): -1: 1
        if func(highs(i)) > func(lows(i))
            if (func(lows(i)) < func(lows(i) - 1e-15)) && (func(highs(i)) > func(highs(i) + 1e-15))
                lows(i) = [];
                highs(i) = [];
            end
        else
            if (func(lows(i)) > func(lows(i) - 1e-15)) && (func(highs(i)) < func(highs(i) + 1e-15))
                lows(i) = [];
                highs(i) = [];
            end
        end
    end
    omega_ans = zeros(1, omega_nums);
    for i = 1: omega_nums
        low = lows(i);
        high = highs(i);
        incline = func(high) > func(low);
        mid = (low + high) / 2;
        limit = 1e-8;
        rounds = 0;
        while (abs(func(mid)) > limit)
            rounds = rounds + 1;
            if rounds > 30
                low = lows(i);
                high = highs(i);
                mid = (low + high) / 2;
                rounds = 0;
                limit = limit * 10;
            end
            if func(mid) > 0
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
        omega_ans(i) = mid;
    end
    %determinant(omega_ans, beta1, theta, H, lambda, P, k1, k2, k3, k4)
end