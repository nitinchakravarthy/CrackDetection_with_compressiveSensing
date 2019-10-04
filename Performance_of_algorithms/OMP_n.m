function [r_data ift_r_data] =  OMP_n(A,b,K)

% Here we recover the compressed data. It should be noted that this
% algorithm should only be executed without clearing any data and should
% immediately follow the compression algorithm
K_target = K ;
opts = [] ;
errorFcn = [] ;
fprintf('Recovery in Progress..... \n');
tic;
[p m] = size(b);
M  = length(A);
r_data = complex(zeros(M,m)) ;
for is = 1:m
    r_data(:,is) = OMP(A, b(:,is), K_target, errorFcn, opts) ;
end
fprintf('Data Recovered \n') ;
t_rec = toc ;
%fprintf('The time taken for the recovery of the data is %d \n', t_rec) ;
ift_r_data = ifft(r_data);

% Now we try to recover the original data in Time Domain
%rec_data = ifft(r_data,'symmetric');
%error = norm(rec_data - data, 2)^2/norm(data,2)^2;
%error = norm(b - A*r_data,1);
%fprintf('The error in the recovery of data(Time Domain) is %e \n', error);

% Now we use the slowMode option in the analysis
% errFcn      = @(a) norm(a-x)/norm(x);
% opts = [];
% opts.printEvery = 25 ;
% opts.slowMode       = true;
% fprintf('Recovery(in slowMode) is in Progress..... \n');
% tic;
% r_data = complex(zeros(N,m)) ;
% for is = 1:m
 %   r_data(:,is) = OMP({Af,At}, b(:,is), K_target, errorFcn, opts) ;
% end
% fprintf('Data Recovered(in slowMode \n') ;
% t_rec1 = toc ;
% fprintf('The time taken for the recovery of the data in slowMode is %d \n', t_rec1) ;

% Now we try to recover the original data in Time Domain
% rec_data = ifft(r_data,'symmetric');
% error1 = norm(rec_data - data, 2)/norm(data,2) ;
% fprintf('The error in the recovery of data(Time Domain)by slowMode is %e \n', error1);

% Note that the code ends by plotting the results.
%figure;
%plot(abs(r_data));