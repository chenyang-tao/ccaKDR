function [B,N] = ccaKDR_CV(X, Y, params)
% ccaKDR     --   ccaKDR with cross validation, Using KNN regression 
%                 to evaluate out-of-sample performance.
%                 Modified from Kenji Fukumizu's gKDR code:
%                 http://www.ism.ac.jp/~fukumizu/gKDR_Matlab.zip
% Input: 
% X         n x p       Predictors
% Y         n x q       Responses
% params    1 x 1       Parameters
% .candx    1 x bx      Kernel bandwidth candidates for X
% .candy    1 x by      Kernel bandwidth candidates for Y]
% .caneps   1 x cx      
% .EK       1 x 1       Reduced dimensionality
% .NCV      1 x 1       
% 
% Output:
% B        p x d       Projection matrix
% N        p x p       Symmetrized gradient space
%
% 
% Version: 1.0
% 
% Copyright (c) Chenyang Tao, 2016. [cytao.fdu(AT)gmail.com]

    candx = params.candx;
    candy = params.candy;
    caneps = params.caneps;
    EK = params.EK;
    NCV = params.NCV;

    N = size(X, 1);

    Method = 'ccaKDR';
    VERBOSE = false;
    
    sgx0=MedianDist(X);   % Basic value for bandwidth for X
    sgy0=MedianDist(Y);   % Basic value for bandwidth for Y
    
    % For cross-validation
%     ridx=randperm(N);  % Random order 
    ridx=[1:N];          % Fixed order here for reproducibility 
    Xr=X(ridx,:);
    Yr=Y(ridx,:);   
    lx=ceil(N/NCV);
    ei=cumsum(lx.*ones(1,NCV),2);
    si=ei-(lx-1).*ones(1,NCV);
    ei(NCV)=N;       % si: staring idx, ei: ending idx
    err_tbl=zeros(length(candx)*length(candy)*length(caneps), NCV);
    
%     fprintf('Cross validation ...');
    
    for h=1:length(candx)
        sgx=sgx0*candx(h);
        for k=1:length(candy)
          sgy=sgy0*candy(k);            
          for ll=1:length(caneps)
            EPS = caneps(ll);
            for i=1:NCV
                ri=si(i):ei(i);
                Xe=Xr; Ye=Yr; 
                Xe(ri,:)=[];
                Ye(ri,:)=[];    % Xe, Ye: trainig sample for CV
                Xt=Xr(ri,:);
                Yt=Yr(ri,:);    % Xt, Yt: test sample for CV
                switch Method
                    case 'ccaKDR'
                        [B t]=ccaKDR(Xe,Ye,EK,sgx,sgy,EPS,EPS);
                    case 'gKDR-i'
                        if M-K <=DEC
                            decd=ones(M-K,1);
                        else
                            dd=floor((M-K)/DEC);
                            r=M-K-dd*DEC;
                            decd=dd*ones(DEC,1);
                            decd(1:r,1)=(dd+1)*ones(r,1);
                        end
                        B=eye(M);
                        for ii=1:length(decd)
                            Ze=Xe*B;
                            Bp=B;
                            dim=M-sum(decd(1:ii,1),1);
                            [B t]=KernelDeriv(Ze,Ye,dim,sgx*sqrt(dim/M),sgy,EPS);
                            B=Bp*B;
                            B=B/sqrtm(B'*B);  
                        end
                    case 'gKDR-v'
                        B=KernelDeriv_var(Xe,Ye,K,sgx,sgy,EPS,50);
                    otherwise
                        error('Error: method mismatch');
                end
                % kNN regression for CV
                B=real(B);
                nnidx=knnsearch(Xe*B,Xt*B, 'K', 5, 'NSMethod', 'kdtree');
    
                Yo=zeros(length(ri),length(Y(1,:)));
                for j=1:length(ri)
                    ii=nnidx(j,:);
                    Yo(j,:)=sum(Ye(ii',:),1)./5;
                end
    
                dd=Yt-Yo;      
                err_tbl((h-1)*length(candy)*length(caneps)+(k-1)*length(caneps)+ll,i)=sum(sum(dd.*dd,1),2)./length(ri);  
        
            end
          end
        end
        if VERBOSE    
            fprintf('.');
            drawnow;
        end
    
    end
    if VERBOSE
        fprintf('\n');
    end
    [c midx]=min(mean(err_tbl,2));
    opth=ceil(midx/(length(candy)*length(caneps)));
    rr=midx-(opth-1)*length(candy)*length(caneps);
    optk=ceil(rr/length(caneps));
    opte=mod(rr,length(caneps));
    if opte==0
        opte=length(caneps);
    end
    
    sgx=sgx0*candx(opth);
    sgy=sgy0*candy(optk);
    EPS=caneps(opte);
    
%     disp([opth optk]);
    
    
    
    
    switch Method
        case 'ccaKDR'
            [B, t, N]=ccaKDR(X,Y,EK,sgx,sgy,EPS,EPS);
%             [B t]=KernelDeriv(X,Y,K,sgx,sgy,EPS);
        case 'gKDR-i'
            error('Not implemented.');
%             B=eye(M);
%             for ii=1:length(decd)
%                 Z=X*B;
%                 Bp=B;
%                 dim=M-sum(decd(1:ii,1),1);
%                 [B t]=KernelDeriv(Z,Y,dim,sgx*sqrt(dim/M),sgy,EPS);
%                 B=Bp*B;
%                 B=B/sqrtm(B'*B);  
%             end
        case 'gKDR-v'     
            error('Not implemented.');
%             B=KernelDeriv_var(X,Y,K,sgx,sgy,EPS,50);
    end
    

end