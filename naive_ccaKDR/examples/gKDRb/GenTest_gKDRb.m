function [X, Y, B0] = GenTest_gKDRb(N)

    Mx = 10;
    My = 1;

    B0 = zeros(Mx, 2);
    B0(1,1) = 1; B0(2,1) = 1;
    B0(1,2) = 1; B0(2,2) = -1;
    B0 = B0/sqrt(2);

    X = 2*rand(N, Mx)-1;
    W = gamrnd(1, 2, [N, My]);
    Z = X*B0;

    Z1 = Z(:,1);
    Z2 = Z(:,2);

    Y = (Z1.^3+Z2).*(Z1-Z2.^3)+W;

end