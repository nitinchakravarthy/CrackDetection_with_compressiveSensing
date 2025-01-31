data = v_acc_end';
L = length(data);
NFFT = 2^nextpow2(L);
ex_data = zeros(NFFT,1);
ex_data(1:L,1) = data(1:L);
Y = fft(ex_data);
Fs = 100;            
T = 1/Fs;             
%mid_length = ceil(L/2);
t = (0:NFFT-1)*T; 
f = Fs*(0:(NFFT/2))/NFFT;
P2 = abs(Y/NFFT);
P1 = P2(1:NFFT/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f_data = Y(1:NFFT/2);
%f_data(1:mid_length) = Y(1:mid_length); 
plot(f,P1)
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%ff_data = f_data;
%n_dash = ceil(length(f_data)/2);
%n_dash = 20000 ;
%ff_data = f_data(1:n_dash,1);
n = length(f_data);
p = 300;
K = 60;
S = 10;
err = 1e-5;
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
ext_x_gomp = zeros(NFFT,1);
ext_x_gomp(1:NFFT/2) = x_gomp(1:NFFT/2);
ext_x_gomp(NFFT/2+1:NFFT) = x_gomp(NFFT/2:-1:1);
ift_x = ifft(ext_x_gomp);
error = norm(data-ift_x(1:L),2)^2/norm(data,2)^2;
fprintf('the error occured in the recovery of the signal:%d\n',error);
fprintf('the time taken by the recovery step using gOMP is:%d\n ',r_time);
