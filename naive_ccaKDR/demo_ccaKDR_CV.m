N = 100;
Mx = 10;
My = 1;

MLOOP = 100;

addpath('examples/gKDRa/');
% addpath('examples/gKDRb/');
% addpath('examples/gKDRc/');
% addpath('examples/gKDRd/');
% Mx = 50;


ccaKDR_CV_opts = struct('candx', [0.25 0.5 0.75 1 2], ...
                        'candy', [0.25 0.5 0.75 1 2], ...
                        'caneps', [1e-5], ...
                        'EK', 1, ...
                        'NCV', 3);

err1 = [];
err2 = [];
te1 = [];
te2 = [];


tid = tic;

for loop=1:MLOOP

    [X, Y, B0] = GenTest_gKDRa(N);
%     [X, Y, B0] = GenTest_gKDRb(N);
%     [X,Y,B0] = GenTest_gKDRc(N,Mx);
%     [X, Y, B0] = GenTest_gKDRd(N);
    
    t0 = toc(tid);
    [B1] = ccaKDR_CV(X, Y, ccaKDR_CV_opts);
    t1 = toc(tid);
    te1(loop) = t1-t0;
    
    err1(loop) = sqrt(trace(B0*B0'*(eye(Mx)-B1*B1'))/trace(B0'*B0));

%     fprintf('ccaKDR:%.3f, gKDR:%.3f\n', err1, err2);
    
    if (mod(loop, 10)==0)
        fprintf('%d/%d\n', loop, MLOOP);
        fprintf('Mean accuracy:%.3f, Time elapsed:%.3fs\n', mean(err1), sum(te1));
    end
    
end

toc(tid);

fprintf('%f, %f, %f\n', mean(err1), std(err1), mean(te1))