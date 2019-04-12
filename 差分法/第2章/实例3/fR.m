function [Rc, Rb] = fR(R_bar, Nc, Nb, K_of_fR, row_vector_fR6)
    cols = size(R_bar, 2);
    R = zeros(2 * Nc + 2 * Nb - 2, cols);
    indexs = [2: Nc, Nc + 2: 2 * Nc, 2 * Nc + 3: 2 * Nc + Nb - 1, 2 * Nc + Nb + 2: 2 * Nc + 2 * Nb - 2];
    indexs_border = [1, Nc + 1, 2 * Nc + 1: 2 * Nc + 2, 2 * Nc + Nb: 2 * Nc + Nb + 1];
    R(indexs, :) = R_bar;
    
    %w
    wcncm1 = R(Nc, :);
    wcncm2 = R(Nc - 1, :);
    wb2 = R(2 * Nc + 3, :);
    wbnbm2 = R(2 * Nc + Nb - 1, :);
    wbnbm3 = R(2 * Nc + Nb - 2, :);
    wbnbm4 = R(2 * Nc + Nb - 3, :);
    swc = sum(R(2: Nc, :));
    
    bw = [zeros(1, cols) ;
        zeros(1, cols);
        zeros(1, cols);
        1/4 * wb2;
        wbnbm3 - 4 * wbnbm2;
        row_vector_fR6 * [swc; wcncm2; wcncm1; wbnbm4; wbnbm3; wbnbm2]];
    
    retw = K_of_fR \ bw;
    R(indexs_border, :) = retw;
    Rc = R(1: 2 * Nc, :);
    Rb = R(2 * Nc + 1: end, :);
end