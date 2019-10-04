% We use the exact ONE-L1 algorithm developed by Zai Yang in the available
% paper. The MATLAB code was taken from his site:
% https://sites.google.com/site/zaiyang0248/publication

% The paper uses the codes to analyse the performance of ONE-L1 algorithms
% on various settings. I am using the eONE-L1 algorithm to compress and
% recover a non-sparse data. I have used FFT to transform the time domain
% signal to frequency domain so as to sparsify it.

[N,m] = size(data);
t = (1 : 1: N) ;
f_data = fft(data);
delta = 0.2 ; % delta = M/N 
rho = 0.2 ; % rho = k/n
M = ceil(delta*N) ;
k = ceil(rho*M) ;
% Now we sparsify the data
pos = zeros(N,m);
for ip = 1: m
  [s, pos(:,ip)] = sort(f_data(:,ip),'descend');
end
s_data = complex(zeros(N,m)) ;
for ip = 1:m
    s_data(pos(1:M),ip) = f_data(pos(1:M),ip);
end

p = randperm(N);
Omega = p(1:M)' ;
% Now we construct the Measurement matrix 'A' and Compress the data
A = pDFT(N,Omega) ;
b = zeros(M,m) ;
fprintf('Compressing the Data......\n');
tic ;
for ic = 1:m
    b(:,ic) = A*s_data(:,ic) ;
end
t_comp = toc ;
fprintf('Done\n');
fprintf('The %d x %d data has been compressed in %d \n', size(data,1), size(data,2),t_comp) ;

r = 1 + delta ;
tol = [1e-5 1e-6] ;
maxiter = [500 500] ;

% Recovery of the data
fprintf('=================== Reconstruction using eONE-L1 ========================\n\n');
r_data = zeros(N,m);
fprintf('Recovery in progress....\n');
tic;
for ir = 1:m
    [r_data(:,ir)] = eONE_L1(A, b(:,ir), r, tol, maxiter, 0);
end
t_rec = toc;
fprintf('Recovery Complete\n');
rec_data = ifft(r_data,'symmetric');
% Now we calculate the error in recovery
error = norm(rec_data - data, 2)/norm(data,2) ;

fprintf('The time taken for recovery of the %d x %d data = %d \n',size(rec_data,1), size(rec_data,2), t_rec);
fprintf('The error in recovery is %e \n',error);

% Now we plot the data and recovered data
figure(1);
subplot(2,1,1);
plot(t,data) ;
title('Original data');
subplot(2,1,2);
plot(t,rec_data);
title('Recovered data');

