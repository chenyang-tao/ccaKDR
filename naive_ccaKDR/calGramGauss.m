function K = calGramGauss(Y, SGY)
    N = size(Y, 1);

    I=eye(N);
    unit=ones(N,N);
    Q=I-unit./N;

    sy2=2*SGY*SGY;
    aa=sum(Y.*Y,2);
    ab=Y*Y';
    D=repmat(aa,1,N);
    yy=max(D + D' - 2*ab, zeros(N,N));
    Gy=exp(-yy./sy2);  
%     Kyo=Q*Gy*Q;
%     Kyo=(Kyo+Kyo')./2;      % Ky(i,j)=Q*exp(-||y(i)-y(j)||^2/SG)*Q
    
    K = Gy;
end