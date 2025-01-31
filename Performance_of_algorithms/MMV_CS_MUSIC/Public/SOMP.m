function [supp,X] = SOMP(Y,A,s,sigma)
% Simulated Orthogona matching pursuit
% works with unknown  s
% by J.C. Ye,  4/5/2012

[m,n] = size(A);
r = size(Y,2);
X= zeros(n,r);



%% noise variance setup
if nargin<4
    sigma =0;  s_max = s;
else
    s_max = m;
end


%% OMP step


supp = [];
supp_c = setdiff(1:n,supp);
res = Y;

while length(supp) < s_max &&  norm(res(:)) > sigma
    supp_c = setdiff(1:n,supp);
    obs = sum(abs(res'*A(:,supp_c)).^2,1);
    [mval,midx] = max(obs);
    supp = union(supp,supp_c(midx));
    [tmpU,R] = qr(A(:,supp),0);
    res= Y - tmpU*tmpU'*Y;
end

supp = sort(supp,'ascend');
X(supp, :) = pinv(A(:,supp))*Y;


