function [supp,r] = MUSIC(Y,A,s,r)
% Standard MUSIC aglorithm
% by J.C. Ye,  4/5/2012

[m,n] = size(A);
N = size(Y,2);
threshold = 0.1;

%% normalize the columns of A
z = sum(abs(A).^2,1).^0.5; 
A = A./z(ones(size(A,1),1),:);


%% subspace estimation from Y*Y'/N
[tV,tD] = eig(Y*Y'/N);
spectrum = diag(tD);
[sval,sidx] = sort(spectrum,'descend');
if nargin < 4
    % estimate r
    if nargin == 3 && size(Y,2) < s
        r = find(sval > 0,1,'last');
    else
        bsval = sval - min(sval);
        tmpv = (bsval(1:m-1)-bsval(2:m))/bsval(1);
        r = find(tmpv > threshold,1,'last');
        if nargin == 3
            r = min([r s]);
        end
    end
    r = min([size(Y,2) r]);
end
U_S = tV(:,sidx(1:r));

%% Standard MUSIC algorithm

supp = [];
supp_c = setdiff(1:n,supp);
obs = sum(abs(U_S'*A(:,supp_c)).^2,1);
[sval,sidx] = sort(obs,'descend');
supp = [supp supp_c(sidx(1:s))];
supp = sort(supp,'ascend');





