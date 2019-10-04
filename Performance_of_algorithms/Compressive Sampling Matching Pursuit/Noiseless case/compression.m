% This code compresses the data imported through the CoSaMP. algorithm.
data = xlsread('data.xls');
data = data(1:10000,:);
[N,m] = size(data);
t = (1: 1: N);
% Now we convert the data to frequency domain
f_data = fft(data);
delta = 0.03; % delta = sampling ratio = M/N
rho = 0.04; % rho = sparsity ratio = K/M
M = round(delta*N) ;
K = round(rho*M) ;
% Now we sparsify the data in frequency domain
s_data = complex(zeros(N,m));
fprintf('Sparsifying the data......\n');
for is = 1:m
    [~,pos] = sort(f_data,'descend');
    s_data(pos(1:K),is) = f_data(pos(1:K),is);
end
fprintf('Done\n');
% Now we try to construct the Measurement data matrix
A = (rand(M,N) + 1i*rand(M,N))/sqrt(2*M) ;
% You can also use the following function instead of the above one
% A = exp(2i*pi*rand(M,N)) ; % Un-comment this and comment the above if you
% wish
% Now we compress the sparsified data
fprintf('Compressing the sparsified data........ \n');
tic;
b = A*s_data ;
t_comp = toc;
fprintf('Done \n')
fprintf('The time taken for compression of the data is %d \n', t_comp);
% This ends the compression of the data.