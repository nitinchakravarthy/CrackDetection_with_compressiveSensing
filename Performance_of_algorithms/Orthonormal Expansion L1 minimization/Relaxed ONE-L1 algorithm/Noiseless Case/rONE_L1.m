function [x_rec, iter] = rONE_L1(A,b,r,tol,maxiter,verbose)

% [x_rec, iter] = rONE_L1(A,b,r,tol,maxiter,verbose)
%
% rONE_L1 solves the ell_1 minimization problem in CS
%   min ||x||_1
%   s.t. Ax = b, where A is partially orthonormal
% Input:
%   A: partially othonormal matrix
%   b: observed sampling data
%   r: exponentially increasing ratio, r > 1
%   tol: convergence tolerance
%   maxiter: maximum iteration
%   verbose: print result if true, and do not if false
% Output:
%   x_rec: reconstructed signal
%   iter: number of iterations
% 
% Written by Zai Yang, 30 Sept 2010
% references: 
% Z. Yang, C. Zhang, J. Deng, and W. Lu, "Orthonornal expansion ...
%     ell_1-minimization algorithms for CS," IEEE Trans. Signal Processing,
%     2011
% Z. Yang, C. Zhang, L. Xie, "On phase transition of CS in the complex
%     domain", IEEE Signal Processing Letters, 2012

if nargin < 6
    verbose = false;
end;
if nargin < 5
    maxiter = 2000;
end;
if nargin < 4
    tol = 1e-5;
end;

z = b;
Atz = A'*z;

n = length(b);
N = length(Atz);

if nargin < 3
    r = min(1 + .04*n/N, 1.02);
end;


x = zeros(N,1);

mu = 1/quantile(abs(Atz),0.99);

k = 0;
converged = false;
while ~converged && k < maxiter
    
    Ax_last = A*x;
    
    x_temp = x + Atz;
    x = sign(x_temp).*max(abs(x_temp) - 1/mu, 0);
    
    z_last = z;
    
    Ax = A*x;
    z = b - (1+1/r)*Ax + Ax_last/r + z_last/r;
    
    Atz = A'*z;
    
    mu = r*mu;
    
    k = k + 1;
    
    converged = (norm(b - Ax)/norm(b) < tol);   % data consistency
    
    if verbose
        fprintf('iter: %d, data consist: %.3e, norm1_x: %.3f\n', k, norm(b - Ax)/norm(b), norm(x,1));
    end;
end;

iter = k;
x_rec = x;