function [ek vk] = fast_eigs(KRA, KRB, KAB, k)
% return the first k eigenvalues of regularized KCCA (KRA, KRB, KAB)
%
% Copyright (c) Chenyang Tao, 2016. [cytao.fdu(AT)gmail.com]
    opts.SYM = true;
    opts.POSDEF = true;
    
    % GPU code
%     gpu_KRA = gpuArray(KRA);
%     gpu_KRB = gpuArray(KRB);
%     gpu_KAB = gpuArray(KAB);
%     gpu_M = gpu_KRA \ gpu_KAB;
%     gpu_M = gpu_KAB' * gpu_M;
%     gpu_M = gpu_KRB \ gpu_M;
%     M = gather(gpu_M);
%     eps = 1e-10;
%     I = eye(size(KRA, 1));
%     KRA = KRA+eps*I;
%     KRB = KRB+eps*I;

    M = linsolve(KRA, KAB, opts);
    M = KAB'*M;
    M = linsolve(KRB, M, opts);
%     M = KRB\(KAB'*(KRA \ KAB));
    
%     tic;
    [v ek] = eigs(M, k);
%     toc;
    ek = diag(ek);
    if (nargout==2)
        vk = v;
    end
end