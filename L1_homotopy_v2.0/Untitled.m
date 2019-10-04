data = v_acc_end';
N = length(data);
M = 150;
T = 60;
delx_mode = 'mil';
A = randn(M,N)/sqrt(M);
f_data = fft(data);
s_data = zeros(size(f_data));
[~,pos] = sort(f_data(:,1),'descend');
s_data(pos(1:T),1) = f_data(pos(1:T),1);
X = cell(2,1);
Y = cell(2,1);
X_out = cell(2,1);

X{1,1} = real(s_data);
X{2,1} = imag(s_data);
Y{1,1} = A*X{1,1};
Y{2,1} = A*X{2,1};
SNR = 100;
%%
in = []; 

    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 1;
%%
for i=1:1:2
    x  = X{i,1};
    y = Y{i,1};
err_fun = @(z) (norm(x-z)/norm(x))^2;
sigma = sqrt(norm(A*x)^2/10^(SNR/10)/M); 
    tau = max(1e-4*max(abs(A'*y)),sigma*sqrt(log(N)));
    tau = tau';
    in.tau = tau;
    in.err_fun = err_fun;
    out = l1homotopy(A,y,in);
    X_out{i,1} = out.x_out;
end
r_s_data =  X_out{1,1} + 1i* X_out{2,1};
r_data = ifft(r_s_data);
error = norm(s_data - r_s_data,2)^2/norm(s_data,2)^2;