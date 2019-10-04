function [x_rec, iter_tot iter] = eONE_L1(A,b,r,Tol,maxIter,verbose)

% [x_rec, iter_tol iter] = eONE_L1(A,b,r,Tol,maxIter,verbose)
%
% eONE_L1 solves the ell_1 minimization problem in CS
%   min ||x||_1
%   s.t. Ax = b,
%  where A is partially orthonormal, see definition in references.
% Input:
%   A: partially othonormal sampling matrix;
%   b: observed sampling data;
%   r: geometricaly increasing ratio, r > 1;
%   Tol: tolerance with 2 entries, 1st refers to the outer loop, 2nd refer
%         to the inner loop;
%   maxIter: maximum iteration with 2 entries, the 1st refers to the outer
%         loop, 2nd refers to the inner loop;
%   verbose: print result if true, and do not if false.
% Output:
%   x_rec: reconstructed signal;
%   iter_tot: total number of (inner) iterations;
%   iter: number of outer iterations.
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
    maxiter = 500;
    maxiterInner = 500;
else
    maxiter = maxIter(1);
    maxiterInner = maxIter(2);
end;
if nargin < 4
    tol = 1e-5;
    tolInner = 1e-6;
else
    tol = Tol(1);
    tolInner = Tol(2);
end;

Atb = A'*b;

N = length(Atb);
n = length(b);

x = zeros(N,1);
y = zeros(n,1);
resid = b;

mu = 1/quantile(abs(Atb),0.99);

iter_tot = 0;
k = 0;
converged = false;
while ~converged && k < maxiter
    
    j = 0;
    convergedInner = false;
    while ~convergedInner && j < maxiterInner
        x_last = x;
        y_temp = resid+y/mu;
        temp_x = x + A'*y_temp;
        x = sign(temp_x).*max(abs(temp_x) - 1/mu, 0);

        resid = b - A*x;
        
        difInner = norm(x-x_last)/(norm(x_last)+1e-10);
        convergedInner = (difInner < tolInner);
        
        j = j + 1;
    end;
    iter_tot = iter_tot + j;
    
    y = mu*y_temp;
    
    mu = r*mu;
    k = k + 1;
    
    converged = (norm(resid)/norm(b) < tol);   % data consistency
    
    if verbose
        fprintf('iter: %d, inner iter: %d, inner consist: %.3e, data consist: %.3e, norm1_x: %.3f\n', k, j, difInner, norm(resid)/norm(b), norm(x,1));
    end;
end;

iter = k;
x_rec = x;