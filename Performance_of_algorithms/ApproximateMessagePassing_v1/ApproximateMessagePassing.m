% Prepare the workspace
%clear; close all; home;

% Dimensions of the problem
data = v_acc_end';
n = length(data);   % Signal length
p = 300;    % Number of measurements
K = 60;     % Number of non-zero elements

% Number of iterations
T = 1000;

% Tolerance
tol = 0.1;

% Generate the problem instance
%A = (1/sqrt(n)).*randn(n, N);
f_data = fft(data);
%n_dash = ceil(length(f_data)/2);
s_data = zeros(length(data),1);
is = 1;
[~,pos] = sort(f_data(:,is),'descend');
s_data(pos(1:K),is) = f_data(pos(1:K),is);
A = (rand(p,n)+1i*rand(p,n))/sqrt(2*p);
x = (s_data);
y = A*x;
% Estimate using Iterative Soft Thresholding
%xist = RecommendedIST(A,y, T, tol);

% Estimate using Approximate Message Passing
tic;
xamp = reconstructAmp(A, y, T, tol, x, 0);
r_time = toc;
% Compute MSE
%erramp = mean((xamp - x).*(xamp - x));
error = norm(xamp-f_data,2)^2/norm(f_data,2)^2;
fprintf('the time taken for recovery of length %d signal is %d\n',n,r_time);
% Print the result
fprintf(' Error in the Recovery using AMP: %.4f\n', error);