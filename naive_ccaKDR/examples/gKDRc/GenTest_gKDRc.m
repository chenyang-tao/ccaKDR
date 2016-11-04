function [X,Y,B0]=GenTest_gKDRc(N,M)
    tid = 1;
    B0 = zeros(M,1);B0(tid)=1;
    AC = 0.5;
    aa=zeros(N,M);
    X=zeros(N,M);
    while nnz(aa)<N*M
        tmpX=randn(N,M);
        X=X.*aa + tmpX.*(1-aa);
        aa=( abs(X) <= 2);
    end
    X=X./2;
    Z=X*B0;    
    Y=(Z(:,1)-AC).^4 .* randn(N,1);
end