function [supp,X] = CSMUSIC(Y,A,s,sigma)
% CSMUSIC   Solve a MMV problem using Compressive MUSIC with  SS-OMP partial support estimation
%
%   CSMUSIC is designed to solve the MMV problem 
%   (MMV)  minimize  ||X||_0  subject to  AX = Y,
%
%   where A is an m-by-n matrix, B is an m-by-N vector with s non-zero rows, and sigma is a
%   noise parameter
%
%   [supp,X] =  CSMUSIC(Y,A,s)       returns estimated support (supp) and
%                                    coefficient X when a sparsity level s is
%                                    known
%   [supp,X] =  CSMUSIC(Y,A,s,sigma) returns estimated support (supp) and
%                                    coefficient X when a sparsity level s is
%                                    unknown but the noise level sigma = norm (Y-AX) is
%                                    known

%   Copyright 2021, Jong Chul Ye, KAIST.
%   http://bisp.kaist.ac.kr/
%  
%   Reference
%   J. M. Kim,O. K. Lee, and J. C. Ye, ?Compressive MUSIC: Revisiting the Link betwee  Compressive Sensing and 
%   Array Signal Processing?, IEEE Trans. on Information Theory, vol. 58, no. 1, pp. 278 - 301, January 2012

[m,n] = size(A);
N = size(Y,2);
X = zeros(n,N);
threshold = 0.1;

% normalize the columns of A
Aorg = A;
z = sum(abs(A).^2,1).^0.5; 
A = A./z(ones(size(A,1),1),:);

% subspace estimation from Y*Y'/N
[tV,tD] = eig(Y*Y'/N);
spectrum = diag(tD);
[sval,sidx] = sort(spectrum,'descend');
bsval = sval - min(sval);
tmpv = (bsval(1:m-1)-bsval(2:m))/bsval(1);
r = find(tmpv > threshold,1,'last');
U_S = tV(:,sidx(1:r));



%% Foward Greey Step
if nargin ==3
    
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
    

    %% Generalized MUSIC step
    [tmpU,R] = qr((eye(m)-U_S*U_S')*A(:,supp),0);
    U_tilde = [U_S tmpU];

    supp_c = setdiff(1:n,supp);
    obs = sum(abs(U_tilde'*A(:,supp_c)).^2,1);
    [sval,sidx] = sort(obs,'descend');
    supp = [supp supp_c(sidx(1:r))];
    supp = sort(supp,'ascend');

    
else
    % unknown sparsity level
    res = U_S;
    supp = [];
    supp_c = setdiff(1:n,supp);
    supp2 = [];
    cag = 1;
    done = 0;
   
    while  ~done && length(supp2) < m
        
        %% Generalized MUSIC step
        U_tilde =orth([U_S A(:,supp)]);
        supp_c = setdiff(1:n,supp);
        obs = sum(abs(U_tilde'*A(:,supp_c)).^2,1);
        [sval,sidx] = sort(obs,'descend');
        supp2 = [supp supp_c(sidx(1:r))];   
        [tmpU,R] = qr(A(:,supp2),0);
        res2= Y - tmpU*tmpU'*Y;

        
        %% one step of OMP
        if norm(res2(:)) > sigma    
            supp_c = setdiff(1:n,supp);     
            [tmpU,R] = qr(A(:,supp),0);
            res= Y - tmpU*tmpU'*Y;
            Phi = res'*A(:,supp_c);
            tmpv = sum(abs(Phi).^2,1);
            [mval,midx] = max(tmpv);
            supp = union(supp,supp_c(midx));
        else
             done = 1;
             supp = supp2;
            
        end
    end

end


supp = sort(supp,'ascend');
X(supp, :) = pinv(Aorg(:,supp))*Y;

