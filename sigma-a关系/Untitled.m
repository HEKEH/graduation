% syms a_sq temp1 temp2 temp3 temp4
% solve (temp1 * a_sq^3 + temp2 * a_sq^2 + temp3 * a_sq - temp4, a_sq)

% -(temp2/(3*temp1))-(2^(1/3)*(-temp2^2+3*temp1*temp3))/(3*temp1(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3)/(32^(1/3)*temp1);
% 
% -(temp2/(3*temp1))+((1+1i*sqrt(3))*(-temp2^2+3*temp1*temp3))/(32^(2/3)*temp1(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))-((1-1i*sqrt(3))*(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))/(62^(1/3)*temp1);
% 
% -(temp2/(3*temp1))+((1-1i*sqrt(3))*(-temp2^2+3*temp1*temp3))/(32^(2/3)*temp1(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))-((1+1i*sqrt(3))*(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4+sqrt(4*(-temp2^2+3*temp1*temp3)^3+(-2*temp2^3+9*temp1*temp2*temp3+27*temp1^2*temp4)^2))^(1/3))/(62^(1/3)*temp1);


sheet1 = [a_cell{p_idxs(1),4001}, a_cell{p_idxs(2),4001}, a_cell{p_idxs(3),4001}];
sheet2 = abs([[phi_c([1/2; 1], p_idxs(1));phi_b([1/2; 1], p_idxs(1))],...
    [phi_c([1/2; 1], p_idxs(2));phi_b([1/2; 1], p_idxs(2))],...
    [phi_c([1/2; 1], p_idxs(3));phi_b([1/2; 1], p_idxs(3))]]);
sheet = zeros(4, 9);
sheet(:, 2: 3: end) = sheet2;
sheet(:, 1: 3: end) = repmat(sheet1, 4, 1);
sheet(:, 3: 3: end) = sheet(:, 2: 3: end) .* sheet(:, 1: 3: end);