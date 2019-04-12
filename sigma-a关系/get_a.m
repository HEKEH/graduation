function a_ans = get_a(p, temp1, temp4, sigma1, omegas, mu_c, mu_b, idx)
    temp2 = 1/4 * p(idx, 2) * sigma1 * omegas(idx);
    temp3 = 1/4 * omegas(idx)^2 * ((p(idx, 1) * mu_c + (1 - p(idx, 1)) * mu_b)^2 + 4 * sigma1^2);
    %方法一
%     syms z
%     a_sqs =  double([root(temp1*z^3 + temp2*z^2 + temp3*z - temp4, z, 1)
%              root(temp1*z^3 + temp2*z^2 + temp3*z - temp4, z, 2)
%              root(temp1*z^3 + temp2*z^2 + temp3*z - temp4, z, 3)]);
%     a_sqs(real(a_sqs) < 0) = [];
%     temp = sqrt(a_sqs);
%     limit = 1e-10;
%     a_ans = temp(abs(imag(temp)) < limit);
%     a_ans = sort(a_ans);

    %方法二
    a_sqs = [-(temp2/(3*temp1))-(2^(1/3)*(-temp2^2+3*temp1*temp3))/(3*temp1*(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3)/(3*2^(1/3)*temp1),...
    -(temp2/(3*temp1))+((1+1i*sqrt(3))*(-temp2^2+3*temp1*temp3))/(3 * 2^(2/3)*temp1*(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))-((1-1i*sqrt(3))*(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))/(6*2^(1/3)*temp1),...
    -(temp2/(3*temp1))+((1-1i*sqrt(3))*(-temp2^2+3*temp1*temp3))/(3 * 2^(2/3)*temp1*(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))-((1+1i*sqrt(3))*(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))/(6*2^(1/3)*temp1)];
    a_sqs(real(a_sqs) < 0) = [];
    temp = sqrt(a_sqs);
    limit = 1e-10;
    a_ans = temp(abs(imag(temp)) < limit);
    count = 0;
    while isempty(a_ans)
        limit = 10 * limit;
        a_ans = temp(abs(imag(temp)) < limit);
        count = count +1;
        if count > 10
            '无解,empty'
            break
        end
    end
    count = 0;
    while length(a_ans) == 2
        limit = limit / 2;
        a_ans = temp(abs(imag(temp)) < limit);
        count = count +1;
        if count > 20
            '无解,两个解'
            break
        end
    end
    a_ans = real(a_ans);
    a_ans = sort(a_ans);


%     da = 5e-3;
%     max_a = 1;
%     a_range = 0: da: max_a; %注意修改
%     a_func = @(a)temp1 * a.^6 + temp2 * a.^4 + temp3 * a.^2 - temp4;
%     ret = a_func(a_range);
%     low_idxs = (ret(1: end - 1) .* ret(2: end)) <= 0;
%     lows = a_range(low_idxs);
%     temp = 1;
%     while length(lows) == 2
%         if temp > 10
%             disp('出错');
%             return;
%         end
%         a_range = 0: da: 2 * max_a;
%         max_a = 2 * max_a;
%         temp = temp + 1;
%         ret = a_func(a_range);
%         low_idxs = (ret(1: end - 1) .* ret(2: end)) <= 0;
%         lows = a_range(low_idxs);
%     end
%     while isempty(lows)
%         if temp > 6
%             disp('出错');
%             return;
%         end
%         a_range = max_a: 10 * da: 10 * max_a;
%         max_a = 10 * max_a;
%         da = 10 * da;
%         temp = temp + 1;
%         ret = a_func(a_range);
%         low_idxs = (ret(1: end - 1) .* ret(2: end)) <= 0;
%         lows = a_range(low_idxs);
%     end
%     highs = lows + da;
%     a_len = length(lows);
%     a_ans = zeros(1, a_len);
%     for i = 1: a_len
%         low = lows(i);
%         high = highs(i);
%         incline = a_func(high) > a_func(low);
%         mid = (low + high) / 2;
%         limit = 1e-20;
%         rounds = 0;
%         while (abs(a_func(mid)) > limit)
%             rounds = rounds + 1;
%             if rounds > 50
%                 low = lows(i);
%                 high = highs(i);
%                 mid = (low + high) / 2;
%                 rounds = 0;
%                 limit = limit * 10;
%             end
%             if a_func(mid) > 0
%                 if incline
%                     high = mid;
%                 else
%                     low = mid;
%                 end
%             else
%                 if incline
%                     low = mid;
%                 else
%                     high = mid;
%                 end
%             end
%             mid = (low + high) / 2;
%         end
%         a_ans(i) = mid;
%     end
end
