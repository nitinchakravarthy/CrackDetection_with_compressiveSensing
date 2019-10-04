function [supp,X] = SeqCSMUSIC(Y,A,s)
% SeqCSMUSIC   Solve a MMV problem using Sequential Compressive MUSIC 
%
%   CSMUSIC is designed to solve the MMV problem 
%   (MMV)  minimize  ||X||_0  subject to  AX = Y,
%
%   where A is an m-by-n matrix, B is an m-by-N vector with s non-zero rows, and sigma^2 is a
%   noise variance
%
%   [supp,X] =  SeqCSMUSIC(Y,A,s)    returns estimated support (supp) and
%                                    coefficient X when a sparsity level s is
%                                    known

%   Copyright 2021, Jong Chul Ye, KAIST.
%   http://bisp.kaist.ac.kr/
% 
%   [1] J. M. Kim, O. K. Lee, and J. C. Ye, Improving Noise Robustness in Subspace-based
%   Joint Sparse Recovery?, IEEE Trans. on Signal Processing , vol. 60, no. 11, pp. 5799-5809, November, 2012


%% initialization
[m,n] = size(A);
N = size(Y,2);
X = zeros(n,N);
threshold = 0.1;

%% normalize the columns of A
Aorg = A;
z = sum(abs(A).^2,1).^0.5; 
A = A./z(ones(size(A,1),1),:);



%% subspace estimation from Y*Y'/N
[tV,tD] = eig(Y*Y'/N);
spectrum = diag(tD);
[sval,sidx] = sort(spectrum,'descend');
bsval = sval - min(sval);
tmpv = (bsval(1:m-1)-bsval(2:m))/bsval(1);
r = find(tmpv > threshold,1,'last');
U_S = tV(:,sidx(1:r));

%%% Initial support estimation routine
    %% partial support estimation using SS-OMP
    supp = [];
    supp_c = setdiff(1:n,supp);
    res = U_S;
    res2 = Y;

    while length(supp) < s-r 
        Phi = res'*A(:,supp_c);
        tmpv = sum(abs(Phi).^2,1);
        [mval,midx] = max(tmpv);
        supp = union(supp,supp_c(midx));
        supp_c(midx) = [];
        tmpU = orth(A(:,supp));
        res = U_S - tmpU*tmpU'*U_S;
        res2 = Y- tmpU*tmpU'*Y;
    end




    %% sequential MUSIC-like step
    while length(supp) < s
        [Utmp,S,V] = svd([U_S A(:,supp)]);
        U_tilde = Utmp(:,1:s);
        supp_c = setdiff(1:n,supp);
        obs = sum(abs(U_tilde'*A(:,supp_c)).^2,1);
        [mval,midx] = max(obs);
        supp = union(supp,supp_c(midx));
        supp_c(midx) = [];      
    end



%% support filtering

supp0 = supp; 
val = zeros(s,1);

for q = 1 : s
    a   = A(:,supp0(q));
    supp_c = setdiff(supp0,supp0(q));
    tmpU   =  orth([U_S,A(:,supp_c)]);
    res = a - tmpU*tmpU'*a;
    val(q) = norm(res);
end
[sval,sidx] = sort(val,'ascend');
supp = supp0(sidx(1:s-r));



%% MUSIC step
while length(supp) < s
        [Utmp,S,V] = svd([U_S A(:,supp)]);
        U_tilde = Utmp(:,1:s);
        supp_c = setdiff(1:n,supp);
        obs = sum(abs(U_tilde'*A(:,supp_c)).^2,1);
        [mval,midx] = max(obs);
        supp = union(supp,supp_c(midx));
        supp_c(midx) = []; 
end

supp = sort(supp,'ascend');
X(supp, :) = pinv(Aorg(:,supp))*Y;


