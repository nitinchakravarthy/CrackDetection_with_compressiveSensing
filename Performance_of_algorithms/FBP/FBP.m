function x = FBP(y, V, alpha, beta, TerMode, eps, K)

% function [x, IndList] = FBP(y, V, K, alpha, beta, TerMode, eps)
%
% INPUT PARAMETERS:
% y:        measurements (Nx1)
% V:        dictionary (MxN)
% alpha:    no. forward selected vectors per iteration
% beta:     no. backwards removed vectors per iteration
% TerMode:  termination mode 
%             ('Err' -> terminate when NMSE<eps)
%             ('K' -> terminate when K vectors selected)
% eps:      error tolerance for termination mode 'Err' 
% K:        sparsity for termination mode 'K'
%
% OUTPUT PARAMETER:
% x:        reconstructed vector

% FBP performs sparse signal reconstruction via the forward-backward pursuit 
% algoritm. 
%
% FBP is an iterative two-stage algorithm. The forward step expands the
% support estimate by alpha atoms, while the backward step removes beta
% atoms from the estimate. Hence, each iteration expands the support
% estimate by alpha-beta. Obviously, alpha>beta must hold.
%
% FBP algorithm:
% initialize support of x: T = 0 
%   iterate until termination:
%     forward step:  add best alpha vectors to T,               
%     backward step: orthogonal projection of y onto T,
%                    remove worst beta vectors from T 
%     update the residue by orthogonal projection of y onto T
%   terminate if 
%       (TerMode = 'Err')   error < eps
%       (TerMode = 'K')     |T| = K
%
% For details, see 
%   - Karahanoglu and Erdogan, "Forward-backward search for compressed
%     sensing signal recovery,", EUSIPCO'2012
%   - Karahanoglu and Erdogan, "Compressed Sensing Signal Recovery via
%     Forward-Backward Pursuit", available as http://arxiv.org/abs/1210.5626, 
%     (submitted to Digital Signal Processing)
%
% Burak Karahanoðlu, Sabanci University, 2012
% email: karahanoglu at sabanciuniv.edu, burak.karahanoglu at gmail.com

Res = y;
M = size(V,1);
N = size(V,2);

if(strcmp(TerMode,'Err'))
    noMaxSupDim = M/2;
    bound = eps*norm(y);
else
    noMaxSupDim = K;  
    bound = 0;
end
T = [];

Vnorm = zeros(N,1);
for k=1:N
    Vnorm(k) = norm(V(:,k));
end;

while norm(Res) > bound % for termination mode 'Err'
   
    % forward step:
    CorrXBasis = (V'*Res)./(Vnorm); % correlation vector
    [a, ind] = sort(abs(CorrXBasis),'descend');
    T = [T; ind(1:alpha)]; % add alpha vectors to support  
    
    % backward step:
    Coefs = V(:,T)\y; % orthogonal projection
    [a, ind] = sort(abs(Coefs));   
    T(ind(1:beta)) = []; % remove the worst beta from support
    
    % update residue by orthogonal projection
    Coefs = V(:,T)\y;    
    Res = y - V(:,T)*Coefs;
    
    if length(T)>=noMaxSupDim % for termination mode 'K'
        break;
    end;
end;

myInd = abs(Coefs)>1e-10; % ignore insignificantly small values
T = T(myInd);
x = zeros(N,1);
x(T) = Coefs(myInd);
