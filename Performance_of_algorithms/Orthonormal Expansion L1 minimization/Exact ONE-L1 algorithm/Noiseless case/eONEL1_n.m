function [x ift_x] = eONEL1_n(A,b,N)
% We recover the data in this file
delta = 0.1;
r = 1 + delta ;
tol = [1e-5 1e-6] ;
maxiter = [500 500] ;

% Recovery of the data
%fprintf('=================== Reconstruction using eONE-L1 ========================\n\n');
[p m] = size(b);
x = zeros(N,m);
%fprintf('Recovery in progress....\n');
for ir = 1:m
    [x(:,ir)] = eONE_L1(A, b(:,ir), r, tol, maxiter, 0);
end
fprintf('Recovery Complete\n');
ift_x = ifft(x,'symmetric');
% Now we calculate the error in recovery
%error = norm(rec_data - data, 2)/norm(data,2) ;

%fprintf('The time taken for recovery of the %d x %d data = %d \n',size(rec_data,1), size(rec_data,2), t_rec);
%fprintf('The error in recovery is %e \n',error);

% NOTE : DoNot clear the data before running the algorithm