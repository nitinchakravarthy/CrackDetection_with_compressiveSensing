function [x ift_x] = test_l1_constraint_n (c_data,n,A) 
%% Basis Pursuit with Douglas Rachford
% Test for DR algorithm for L1 minimization (BP).
% We do here a compressed sensing resolution
% (random matrix).

%%
% Add the toolbox.

%addpath('../');
%addpath('../toolbox/');

%% 
% Dimensionality of the signal and number of measurements.
y = c_data;
[p m] = size(y);
%A = (rand(p,n) + 1i* rand(p,n))/sqrt(2*p);
%%
% We aim at solving 

%%
% |min_{A*x=y} norm(x,1)|

%%
% This can be rewriten as the minimization of |F(x)+G(x)|
% where |F=norm(x,1)| and |G=i_{A*x=y}| is the indicator function.


%%
% The proximity operator of the L1 norm is the soft thresholding.

ProxF = @(x,tau)perform_soft_thresholding(x, tau);

%%
% The proximity operator of the indicator of |A*x=y| is the orthogonal
% projection on A*x=y.

pA = A'*(A*A')^(-1);
ProxG = @(x,tau)x + pA*(y-A*x);

%%
% Create a function to record the values of F and the constraint at each iteration.

F = @(x)norm(x,1);
Constr = @(x)1/2*norm(y-A*x)^2;
options.report = @(x)struct('F', F(x), 'Constr', Constr(x));

%%
% Run the algorithm. 

options.gamma = 5;
options.niter = 5000;
[x,R] = perform_dr(zeros(n,m), ProxF, ProxG, options);

%%
% Display the solution. At convergence, it should be of sparsity |p|.
%{
figure(2);
plot(x);
axis tight;
%}
%%
% Retrieve the F and constraint function values.

f = s2v(R,'F');
constr = s2v(R,'Constr');

%%
% Display.

%{
figure(3);
subplot(2,1,1);
plot(f(2:end));
axis tight; title('Objective');
subplot(2,1,2);
plot(constr(2:end));
axis tight; title('Constraint');
rec_data = ifft(x);
%}
%{
figure;
plot(1:length(x),abs(x))
title('Recovered Sparse Data')
xlabel('------Time-----' )
ylabel('----Sparse Data----');
%}
ift_x = ifft(x);

%figure(2);
%plot(1:length(data),real(r));
%{
d = data;
error = norm(d-r,2)^2/norm(d,2)^2;
fprintf('the error obtained in the recovery is%d\n', error );
fprintf('time taken is :%d\n',time);
error_com = norm(d -ifft(f_data),2)^2/norm(d,2)^2;
fprintf('the error obtained in the compression is%d\n', error_com );
%}