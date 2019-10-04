%basic algorithm 
%   ROBUST PRINCIPAL COMPONENTS ANALYSIS BASED ON PROJECTION PURSUIT
%
% Croux and Ruiz-Gazen
% High Breakdown estimators for Principal Components: the
% Projection-Pursuit Approach Revisited, Journal of Multivariate Analysis,
% forthcoming.
%
% 
% {lambda,eigv,scores}=robpca(x,pp,s)
%
%     x: data-matrix  (n,p)
%     pp: number of desired eigenvectors  (pp<=p)
%     s: name of scale estimator, e.g. mad
%     lambda: eigenvalues
%     eigv:  eigenvectors
%     scores: matrix of scores (n,pp)
%
%
% uses: L1median.m  (to compute the spatial median)
%       s.m  (procedure to compute the robust scale estimator, e.g. mad.m)
%
%%%%%%%%%%%%%%%%%%%%%%
%Example, where x is a datamatrix with 500 observations and 20 variables:
%	  x=randn(500,20);
%	  [lambda1,eigv1,scores1]=robpca(x,2,@mad);
%	  [lambda2,eigv2,scores2]=robpca(x,2,@qn);
%
% As such, lambda1 contains the "eigenvalues" using the projection-pursuit 
% approach with the Median Absolute Deviation scale estimator, using the Croux 
% and Ruiz-Gazen algorithm. The eigenvalues are in eigv1, and the scores
% (projections of the data on the principal components) in scores1.
%%%%%%%%%%%%%%%%%%%%%%

function [lambda,veig,scores]=robpca(x,pp,s);


[n,p]=size(x);

if (pp > min([n p])); error('pp too large'); end;

if (p>n);
[v,d,u]=svd(x',0);
x=u*d;
pold=p;
p=n;
else
pold=p;
end;

m=L1median(x);
y=x-repmat(m,n,1);
bigscores=zeros(n,pp);

veig=[];lambda=[];

for k=1:pp;
  if (k < p);
   pcol=zeros(n,1);
   for i=1:n;
     pyi=y(i,:); pyi=pyi';
     npyi=norm(pyi);
     if (npyi==0);
     pcol(i)=0;
     else
     pyi=pyi/npyi;
     pcol(i)=feval(s,y*pyi);
     end
   end
   [lambdastar,istar]=max(pcol);
   lambda=[lambda ; lambdastar];
   vhelp=y(istar,:)';
   vhelp=vhelp/norm(vhelp);
   scores=y*vhelp;
   y=y-scores*(vhelp)';
  else % last eigenvector is automatically found
   i=1;
   while (norm(y(i,:))==0),
   i=i+1;
   end
   vhelp=y(i,:)';
   vhelp=vhelp/norm(vhelp);
   scores=y*vhelp;
   lambda=[lambda ; feval(s,y*vhelp)];
  end
 veig=[veig vhelp];
 bigscores(:,k)=scores;
k=k+1;
end;

if (pold>n);
veig=v*veig;
end;

lambda=lambda.^2;
scores=bigscores;



