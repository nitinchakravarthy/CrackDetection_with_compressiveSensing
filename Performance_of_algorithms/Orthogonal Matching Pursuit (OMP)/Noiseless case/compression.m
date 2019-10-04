% Here we compress the data by transforming it to Frequency basis (using
% Fourier transform).
% Case : NOISELESS
%data = xlsread('data.xls');
data = v_acc_end;
[N,m] = size(data);
t = (1 : 1: N) ;
f_data = fft(data);
delta = 0; % delta = M/N i.e. the sampling ratio.
rho = 0; % rho is the sparsity ratio i.e. rho = K/M
M = 120 ;
K = 20 ;
s_data = complex(zeros(N,m));
for is = 1:m
    [~,pos] = sort(f_data(:,is),'descend');
    s_data(pos(1:K),is) = f_data(pos(1:K),is);
end
% Now we construct the measurement matrix 'A'
A =  (randn(M,N) + 1i*randn(M,N))/sqrt(2*M);
% It is not necessary that the measurement matrix 'A' has unit norm columns
% but it helps if it does.
Af = @(x) A*x ;
At = @(x) A'*x ;
% Now we compress the data.
fprintf('Compressing the data.......\n');
tic ;
b = complex(zeros(M,m));
for is = 1:m
    b(:,is) = A*s_data(:,is) ;
end
fprintf('Compression Complete \n');
t_comp = toc;
fprintf('The time for compression is %d \n',t_comp);