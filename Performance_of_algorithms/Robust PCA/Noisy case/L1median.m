function [mX]=L1median(X,tol);

% L1MEDIAN calculates the multivariate L1 median 
% I/O: [mX]=L1median(X,tol);
%
% X is the data matrix 
% tol is the convergence criterium; the iterative proces stops when ||m_k - m_{k+1}|| < tol.
%
% Ref: Hossjer and Croux (1995) "Generalizing Univariate Signed Rank Statistics for Testing
% and Estimating a Multivariate Location Parameter", Non-parametric Statistics, 4, 293-308.
% Translated from the Gauss code of Hossjer and Croux (1995) in matlab by Sabine Verboven, Antwerp University. 

if nargin <2
   tol=1.e-08;
end;
[n,p]=size(X);
maxstep=200;
%initializing starting value for m
m=median(X);
k=1;
while (k<=maxstep)
   mold=m;
   Xext=sortrows([norme(X-repmat(m,n,1)) X],1);
   dx=Xext(:,1);
   X=Xext(:,2:p+1);
   if all(dx)
      w=1./dx;
   else
      ww=dx(all(dx,2));
      w=1./ww;
      w=[zeros(length(dx)-length(w),1);w];
   end
   delta=sum((X-repmat(m,n,1)).*repmat(w,1,p),1)./sum(w);
   nd=norme(delta);
   if all(nd<tol)
      maxhalf=0;
   else
      maxhalf=log2(nd/tol);
   end
	m=mold+delta;   %computation of a new estimate
   nstep=0;
   while all(mrobj(X,m)>=mrobj(X,mold))&&(nstep<=maxhalf)
      nstep=nstep+1;
      m=mold+delta./(2^nstep);
   end
   if (nstep>maxhalf)
      mX=mold;
      break
   end
   k=k+1;
end
if k>maxstep
   display('Iteration failed')
end
mX=m;



%-----
function n=norme(X)

%NORME calculates the euclidian norm of matrix X
% the output is a columnvector containing the norm of each row
%I/O: n=norme(X);

n = sqrt(sum(X.^2,2));

%--------
function s=mrobj(X,m)

%MROBJ computes objective function in m based on X and a

xm=norme(X-repmat(m,size(X,1),1));
s=sum(xm,1)';

