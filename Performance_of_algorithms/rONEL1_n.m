function [x ift_x] = rONEL1_n(A,b,N)
% Now we recover the data using the relaxed ONE-L1 algorithm
delta = .1;
r = min(1+0.04*delta, 1.02) ;
tol = 1e-5 ;
maxiter = 1000 ; % The maximum number of iteration to be done
%fprintf('=================== Reconstruction using rONE-L1 ========================\n');
%fprintf('Recovering the Data............\n');
[p m] = size(b);
x = zeros(N,m);
for i = 1:m
    [x(:,i)] = rONE_L1(A, b(:,i), r, tol, maxiter, 0);
end
fprintf('Done\n');
ift_x = ifft(x, 'symmetric');
%error = norm(rec_data - data,2)/norm(data,2) ;
%fprintf('Relative root mean squared error: %e \n', error);
%fprintf('The time taken for compression of data is %d seconds\n', r_time);