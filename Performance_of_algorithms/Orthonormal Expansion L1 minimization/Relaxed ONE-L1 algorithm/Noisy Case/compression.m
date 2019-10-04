% We compress the data here.
% We attempt to use relaxed ONE-L1 algorithm to reconstruct the data. The
% code was based on the methods discussed in the following paper

% Orthonormal Expansion L1-Minimization Algorithms for Compressed Sensing %
%         By: Zai Yang, Cishen Zhang, Jun Deng, and Wenmiao Lu

% The original codes was for sparse signals but I have modified it so as to
% use for Non-sparse data.
%data = xlsread('data.xls');
data = v_acc_end';
[N,m] = size(data);
noise = randn(N,m);
noisy = data  ;
t = (1 : 1: N);
new_data = fft(noisy);
% Now we sparsify the above fft data so as to compress it.
delta = 0.01 ; % Change the value for changing the compression ratios
rho = 0.2 ;

n = ceil(delta*N) ;
K = ceil(rho*n) ;


% Sparsifying. The process is not perfect as I am only taking the highest
% 'rho*n' elements. You can also try to isolate the local maximas and then
% take the mean and standard deviation of the recovered data so as to find
% a Threshold and then use it to zero the data.
s_data = complex(zeros(N,m)) ;
for i = 1:m 
    [~, pos] = sort(new_data(:,i),'descend');
    s_data(pos(1:K),i) = new_data(pos(1:K),i) ;
end
    

% Now we create the sensing matrix 'A' using the pDFT(partial Fourier
% Transformation) function. The function is also developed by Zai Yang and
% the sensing matrix created using this function has its own class i.e.
% pDFT.
p = randperm(N);
Omega = p(1:n)';

A = pDFT(N,Omega);

% Now we compress each and every column of the sparsified matrix
b = zeros(n,m) ;
fprintf('Compressing the Data........\n');
s_data = new_data;
tic ;
for i = 1:m
    b(:,i) = A*s_data(:,i);
end
c_time = toc ;
fprintf('Done\n');
fprintf('The time taken for compression of %dx%d data is %d \n\n', size(s_data,1),size(s_data,2),c_time) ;
