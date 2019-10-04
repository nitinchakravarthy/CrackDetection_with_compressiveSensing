function [r_data ift_r_data] = CoSaMP_n(A,b,K)
% Now we try to recover the sparsified data using the CoSaMP algorithm
%{
Comments:
    This is sensitive to K_target
    If the signal is noiseless, then K_target should be an overestimate
    of the true sparsity level.
    But if the signal is noisy, then the optimal performance usually
    happens if K_target is an underestimate of the true sparsity level.
    (This rule-of-thumb is also true for OMP )

    Also, convergence may be very slow if K_target is large.
    Furthermore, each iteration is slower if K_target is large.

    The "HSS" mode may help for noisy vectors.

    So, summary of recommendations for the case with noisy data:
    "HSS" and "two_solves" mode should be "true"
    "K_target" should be small, like 5% of N
    "addK" should be like 1*K_target or less, not 2*K_target.

%}
opts            = [];
opts.maxiter    = 50;
opts.tol        = 1e-8;
opts.HSS        = true;
opts.two_solves = true; % this can help, but no longer always works "perfectly" on noiseless data
opts.printEvery = 10;
errFcn = [] ;
%tic;
fprintf('Recovery in progress........ \n')
N= length(A);
[p m] = size(b);
r_data = complex(zeros(N,m));
for is = 1:m
    [r_data(:,is)] = CoSaMP(A, b(:,is), K,errFcn, opts);
end
ift_r_data = ifft(r_data);
%t_rec = toc;
%fprintf('Recovery Complete\n');
%fprintf('The time taken for the recovery of the data is %d \n', t_rec);

% Now we try to recover the original data using inverse fft
%rec_data = ifft(r_data,'symmetric');
%error = norm(rec_data-data,2)^2/norm(data,2)^2 ;
%fprintf('The error in recovery is %e \n\n', error);
