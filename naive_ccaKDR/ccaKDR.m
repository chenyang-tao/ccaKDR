function [Bx,t,Nx] = ccaKDR(X,Y,d,SGX,SGY,EPSX,EPSY)
% ccaKDR     --   Canonical kernel dimension reduction
%                 Naive implementation
%
% Input: 
% X         n x p       Predictors
% Y         n x q       Responses
% d         1 x 1       Reduced dimensionality
% SGX,SGY   1 x 1       Kernel bandwidth parameter
% EPSX,EPSY 1 x 1       Regularization parameter
% 
% Output:
% Bx        p x d       Projection matrix
% t         d x 1       Sum of d-leading canonical correlations
% Nx        p x p       Symmetrized gradient space
%
% 
% Version: 1.0
% 
% Copyright (c) Chenyang Tao, 2016. [cytao.fdu(AT)gmail.com]
%
% REF:  Tao C. and Feng J. (2016) Canonical kernel dimension reduction


    [N,Mx] = size(X);
    
    Q=eye(N)-1/N*ones(N);
    Gx=calGramGauss(X, SGX);
    Gy=calGramGauss(Y, SGY);
    Kx=Q'*Gx*Q;
    Ky=Q'*Gy*Q;

    % Symmetrize
    Kx=(Kx+Kx')/2;
    Ky=(Ky+Ky')/2;

    k = d;

    [Wx, Wy, EV] = getKCCA2(Kx, Ky, EPSX, EPSY, k);

    Dx=reshape(repmat(X,N,1),N,N,Mx);
    Xij=Dx-permute(Dx,[2 1 3]);
    Xij=Xij./SGX/SGX;
    Hx=Xij.*repmat(Gx,[1 1 Mx]);
%     Hx=Xij.*repmat(Kx,[1 1 Mx]);
    Hx=permute(Hx,[1,3,2]);
    Hx=reshape(Hx,N*Mx,N);
    Xd=Hx*Wx;       % N*Mx \times k
    EV=repmat(EV',N*Mx,1);
%     EV=repmat(ones(1,k),N*Mx,1);
    Xd=Xd.*EV;      % N*Mx \times k
    Xd=reshape(Xd,N,Mx,k);
    Xd=permute(Xd,[1,3,2]); % N*k \times Mx
    Xd=reshape(Xd,N*k,Mx);
    
    % Solution
    Nx=Xd'*Xd;
    [Vx,Ex]=eig(Nx);
    [e,idx]=sort(diag(Ex),'descend');
    t=sum(e(idx(1:d)));

    Bx=Vx(:,idx(1:d));
    
end