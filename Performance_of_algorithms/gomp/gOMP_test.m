data = Data(1:10000,3);
[n m] = size(data);
p = 2000;
K = 1000;
err    = 1e-15;
S = 1;
f_data = fft(data);
%{
ff_data = f_data;
%n_dash = ceil(length(f_data)/2);
%n_dash = 20000 ;
%ff_data = f_data(1:n_dash,1);
s_data = zeros(length(f_data),1);

is = 1;
[~,pos] = sort(f_data(:,is),'descend');
s_data(pos(1:K),is) = f_data(pos(1:K),is);
A = (rand(p,n)+1i*rand(p,n))/sqrt(2*p);
x = (s_data);
y = A*x;
tic;
[x_gomp, support, iteration] = gomp(y, A, K, S, err);
r_time = toc;
error = norm(data-ifft(x_gomp),2)^2/norm(data,2)^2;
fprintf('the error occured in the recovery of the signal:%d\n',error);
fprintf('the time taken by the recovery step using gOMP is:%d\n ',r_time);
%}