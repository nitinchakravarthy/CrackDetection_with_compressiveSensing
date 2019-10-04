% Now we recover the data using the relaxed ONE-L1 algorithm
r = min(1+0.04*delta, 1.02) ;
tol = 1e-5 ;
maxiter = 1000 ; % The maximum number of iteration to be done
fprintf('=================== Reconstruction using rONE-L1 ========================\n');
fprintf('Recovering the Data............\n');
tic ;
r_data = zeros(N,m) ;
for i = 1:m
    [r_data(:,i)] = rONE_L1(A, b(:,i), r, tol, maxiter, 0);
end
r_time = toc ;
fprintf('Done\n');
rec_data = ifft(r_data, 'symmetric');
error = norm(rec_data - data,2)/norm(data,2) ;
fprintf('Relative root mean squared error: %e \n', error);
fprintf('The time taken for compression of data is %d seconds\n', r_time);