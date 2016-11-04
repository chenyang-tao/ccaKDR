function [X, Y, B0]=GenTest_gKDRd(N)
    % N=1000,2000
    Mx = 50;
    My = 1;
    d=1;

    B0= zeros(Mx, 2*d);
    for i=1:d
        B0(2*i-1:2*i,2*i-1:2*i) = [1, 1; 1, -1]/sqrt(2);
    end

    X = 2*rand(N, Mx)-1;
    W = laprnd(N, My, 0, 2);
    Z = X*B0;

    Y = zeros(N, My);
    for i=1:d
        Z1 = X*B0(:, 2*i-1);
        Z2 = X*B0(:, 2*i);
        Y = Y+(Z1.^3+Z2).*(Z1-Z2.^3);
    end
    Y = Y+W;

end