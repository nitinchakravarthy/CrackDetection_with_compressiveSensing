function [ Y,A,X, supp] = forward_gaussian( m, n, s, r, snr, N, a, tau )


%% Gaussian sensing matrix
    A = randn(m, n) *sqrt(1/m)+a;
    z = sum(abs(A).^2,1).^0.5; 
    A = A./z(ones(size(A,1),1),:);
    
%% Noise generation

    flag_complex = 0;
    
    rp = randperm(n);
    supp = rp(1:s);
    supp = sort(supp,'ascend');

    % mixed multichannel model
    Psi = zeros(n,s);
    for k = 1:s
        Psi(supp(k),k) = 1;
    end
    
    
    if flag_complex == 1
        U = randn(s,r) + j*randn(s,r);
    else
        U = randn(s,r);
    end
    U = orth(U);
    Lambda = diag(tau.^(0:r-1));%eye(r);
    if flag_complex == 1
        V = randn(N,r)/sqrt(2*N) + j*randn(N,r)/sqrt(2*N);
    else
        V = randn(N,r)/sqrt(N);
    end    
    X0 = Psi*U*Lambda*V';
    Y = A*X0;
    
    if flag_complex == 1
        Z = randn(m,N) + j*randn(m,N);
    else
        Z = randn(m,N);
    end        
    Z = Z / norm(Z, 'fro') * norm(Y,'fro') * 10^(-snr/20);
    Y = Y + Z;
%     if ~isinf(snr)
%         Y = Y + Z;
%     end



% %%  Signal generation
%     tmpv = randperm(n);
%     supp  = tmpv(1:k);
%     
%     Psi = orth(rand(k,r));
%     Phi = orth(rand(N,r));
%    
%     X = zeros(n, N);
%     X(supp, :) = Psi*Phi';
%     
%  %% Data
%     Y = A*X;            
%     val = sqrt((Y(:)'*Y(:))* 10^(-snr/10));
%     W = randn(size(Y));
%     tmp = sqrt(W(:)'*W(:));
%     Y = Y + W/tmp*val;

X = X0;
        
end

