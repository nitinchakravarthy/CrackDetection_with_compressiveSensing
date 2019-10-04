function [x_gomp, support, iteration] = gomp(y, A, K, S, err)

% Generalized Orthogonal Matching Pursuit (gOMP) is a greedy algorothm that
% provides approximate solution to the problem: min ||x||_0 such that Ax = y. 
% gOMP extends the conventional Orthogonal Matching Pursuit (OMP) algorithm
% by allowing S indices to be chosen per iteration.   

% Input:
% y		    : measurements
% A       	: measurement matrix
% K	    	: Sparsity level of underlying signal to be recovered
% S	    	: number of indices chosen per iteration
% err       : residual tolerant


% Output:
% x_gomp    : estimated sparse signal
% support   : estimated support set
% iteration : # of performed iterations 

% Copyright (c) Jian Wang, Seokbeop Kwon and Byonghyo Shim, 2012

 if  nargin < 5
	   err    = 1e-5;
 end 
    
	x_gomp	  = zeros(size(A,2), 1);
	residual  = y;
	supp	  = [];
	iteration = 0; 
	
	while (norm(residual) > err && iteration < min(K, floor(size(A,1)/S))) 
		   iteration          = iteration + 1;
		   [~, idx]           = sort(abs(A' * residual), 'descend');
		   supp_temp          = union(supp, idx(1:S));

	   if (length(supp_temp) ~= length(supp))
           supp	              = supp_temp;
		   x_hat			  = A(:,supp)\y;
		   residual           = y - A(:,supp) * x_hat; 
       else
		   break;
       end
    end
    
 	x_gomp(supp)	          = A(:,supp)\y;
	[~, supp_idx]             = sort(abs(x_gomp), 'descend');
	support                   = supp_idx(1:K); 
	x_gomp                    = zeros(size(A,2), 1);
    x_gomp(support)           = A(:,support)\y;
end
