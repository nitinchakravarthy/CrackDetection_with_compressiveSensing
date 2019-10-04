% We compress the data here using the exact Orthonormal Expansion L1
% minimization algorithm
% We use the exact ONE-L1 algorithm developed by Zai Yang in the available
% paper. The MATLAB code was taken from his site:
% https://sites.google.com/site/zaiyang0248/publication

% The paper uses the codes to analyse the performance of ONE-L1 algorithms
% on various settings. I am using the eONE-L1 algorithm to compress and
% recover a non-sparse data. I have used FFT to transform the time domain
% signal to frequency domain so as to sparsify it.
data = xlsread('data.xls');
data = data(1:10000,:);
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