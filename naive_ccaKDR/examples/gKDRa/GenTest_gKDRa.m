function [X, Y, B0] = GenTest_gKDRa(N)

    Mx = 10;
    My = 1;

    SIG = 1e-1;

    B0 = zeros(Mx, 1);
    B0(1) = 1;
    B0(2) = 2;
    B0 = B0/sqrt(5);

    X = 2*rand(N, Mx)-1;
    W = SIG*randn(N, My);
    Z = X*B0;

    Y = Z.*sin(Z)+W;
    
end