%%
data = Data(1:10000,2:1:17);
K= 400;
p = 800;
[n m] = size(data); 
f_data = fft(data);
%%
s_data = zeros(n,m);

[~,pos] = sort(f_data(:,3),'descend');
s_data(pos(1:K),3) = f_data(pos(1:K),3);
 for is = 1:m
 s_data(pos(1:K),is) = f_data(pos(1:K),is);
 end
 %%
X = s_data;
A = (randn(p,n)+ 1i*randn(p,n))/sqrt(2*p);
Y = A*X;
r_time = zeros(4,1);
error = zeros(4,1);
Xh = cell(4,1);
%%
tic;
[actSet, Xh{1,1}] = SeqCSMUSIC(Y, A, K);
r_time(1) = toc;
 error(1) = norm(data-ifft(Xh),2)/norm(data,2);
 %%
 tic;
[actSet, Xh{2,1}] = CSMUSIC(Y, A, K);
r_time(1) = toc;
 error(1) = norm(data-ifft(Xh),2)/norm(data,2);
 %%
 tic;
[actSet, Xh{3,1}] = SOMP(Y, A, K);
r_time(1) = toc;
 error(1) = norm(data-ifft(Xh),2)/norm(data,2);
 tic;
 %%
[actSet, Xh{4,1}] = MUSIC(Y, A, K);
r_time(1) = toc;
 error(1) = norm(data-ifft(Xh),2)/norm(data,2);
 %%