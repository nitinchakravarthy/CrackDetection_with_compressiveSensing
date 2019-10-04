
function [c_data, K, f_data, s_data,ff_data, A,p] = compressor(data,method_comp,comp_type,noise_ind)
%data = data(:,:);
%rho = 0.1;   % p/n = rho
%delta = 0.3; % k/p  = delta
%p = rho * n;
%K = delta * p;
%%
if comp_type == 4
p = 800;
else
p = 1200;
end
K =400;
 m = size(data,2);
% Sensing matrix.
ff_data = fft(data);
f_data = ff_data(1:(length(ff_data)/2)+1,1:m);
n  = length(f_data);
s_data = zeros(n,m);
%%

if comp_type 
[~,pos] = sort(f_data(:,3),'descend');
s_data(pos(1:K),3) = f_data(pos(1:K),3);
 for is = 1:m
 s_data(pos(1:K),is) = f_data(pos(1:K),is);
 end
 
else
for is = 1:m
[~,pos] = sort(f_data(:,is),'descend');
s_data(pos(1:K),is) = f_data(pos(1:K),is);
end

end
%%
if (~(method_comp == 5 || method_comp == 6 )) %%%%% For eONEL1 and rONEL1 Projection Mtraix has to be Orthonormal
    if method_comp == 8
    A = randn(p,n)/sqrt(p);  %%%% AMP requires real valued projection matrix
    else
    A = (rand(p,n)+1i*rand(p,n))/sqrt(2*p);
    end
    c_data = A*s_data;
    c_noisy_data = awgn(c_data,10);
elseif(method_comp == 5 || method_comp == 6)%%%%% For eONEL1 and rONEL1 Projection Mtraix has to be Orthonormal
        perm = randperm(n);
        Omega = perm(1:p)' ;
    % Now we construct the Measurement matrix 'A' and Compress the data
    A = pDFT(n,Omega) ;
    b = zeros(p,m) ;
    fprintf('Compressing the Data......\n');
    tic ;
    for ic = 1:m
        b(:,ic) = A*s_data(:,ic) ;
    end
    c_data = b;
    c_noisy_data = awgn(c_data,10);
if noise_ind
   c_data = c_noisy_data;
end
end
