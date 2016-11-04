function [W1, W2, EV] = getKCCA2(K1, K2, EPS1, EPS2, k)
% Compute the solution of a regularized KCCA problem
%
% Copyright (c) Chenyang Tao, 2016. [cytao.fdu(AT)gmail.com]
    n = size(K1, 1);
    I = eye(n);

    RK1 = (K1 + n*EPS1*I)^2;
    RK2 = (K2 + n*EPS2*I)^2;
    K12 = K1*K2;

    [ek vk] = fast_eigs(RK1, RK2, K12, k);
    W2 = vk;
    W1 = RK1 \ (K12*W2);
    if (nargout==3)
        EV = ek;
    end
end