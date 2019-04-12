function [K_of_fR, row_vector_fR6] = get_params_of_fR(beta, lambda, H, theta, m, dxc, dxb)
    k1 =  m / (beta * cos(theta));
    k2 =  m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta * cos(theta));
    k3 =  m * (- 2 * H * tan(theta) + 1)^2 * lambda^2 * cos(theta) / (4 * beta);
    k4 =  -m / (beta * cos(theta));
    
    l1 = k2 * dxb^3 * dxc;
    l2 = k2*dxb^3*dxc + (3*k1*dxb^3)/dxc;
    l3 = 18 - 4 * k4 * dxb^2;
    l4 = dxb^2 * (2 * k3 * dxb + 3 * k4) - 5;
    r1 = -2 * k2 * dxb^3 * dxc;
    r2 = -((k1 * dxb^3)/dxc);
    r3 = (4 * k1 * dxb^3) / dxc;
    r4 = 3;
    r5 = -14;
    r6 = 24 - k4 * dxb^2;
    
    K_of_fR = ...
    [1,     0,     0,     0,     0,     0;
     0,     1,     0,     0,     0,     -cos(theta)^2;
     0,     0,     1,     0,     0,     0;
     0,     0,     -3/4,     1,     0,     0;
     0,     0,     0,     0,     -5,     2;
     l1,    l2,     0,     0,     l3,    l4;];
    
    row_vector_fR6 = [r1, r2, r3, r4, r5, r6];
end