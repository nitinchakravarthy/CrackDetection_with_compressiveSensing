% script to demonstrate use of FBP function.
%  - creates a random sparse vector x and a random observation matrix Phi.
%  - computes random observation of x as y = Phi*x
%  - reconstructs x from y via forward-backward pursuit (FBP).
%
% Burak Karahanoðlu, Sabanci University, 2012
% email: karahanoglu at sabanciuniv.edu, burak.karahanoglu at gmail.com

clc
%% create test data 
% data properties
data = v_acc_end';
x= fft(data);
N = length(data);
M = 400;
K = 60;
%{
K = 30;     % sparsity
N = 256;    % # all dimensions 
M = 100;    % # observations

% create a random sparse vector
display('Creating random sparse test vector x...');
x = zeros(N,1);
x_ind = randperm(N);
x_ind = x_ind(1:K);
x(x_ind) = randn(K,1);          % standard Gaussian entries
%}
% create random observation matrix
display('Creating random observation matrix Phi...');
%Phi = randn(M,N)/N;
% % you may want to normalize Phi! (usually improves recovery with FBP)
% for k=1:N
%     Phi(:,k) = Phi(:,k)/norm(Phi(:,k)); % normalize Phi
% end;
% % or,
% Phi = orth(Phi')';
Phi = (randn(M,N) + 1i*randn(M,N))/sqrt(2*M);

y = Phi*x; % observed vector

%% run FBP
eps = 1e-7;
fs = 10;
br = 7;
TerMode = 'Err';
display('Running FBP...');
xhat = FBP(y, Phi, fs, br, TerMode, eps, K);

%% display results
%{
xhat_ind = find(xhat);
display('RECONSTRUCTION RESULT:')
x_ind = sort(x_ind,'ascend')';
if (isequal(x_ind, xhat_ind))
    display(' x is exactly reconstructed.');
else
    noIdentifiedComp = sum(ismember(xhat_ind, x_ind));
    display(' Exact reconstruction failed:') 
    display(['   ' num2str(noIdentifiedComp) ' entires of x were correctly reconstructed.']); 
    display(['   ' num2str(K-noIdentifiedComp) ' entires of x were missed.']); 
end
%}
%% plot signals 
%MSE = mean((x-xhat).^2);
%display([' Mean Squared Error: ' num2str(MSE)]);
%figure; plot(1:N,x,'--b'); hold on; plot(1:N,xhat,'-.r'); plot(1:N,(x-xhat),':g');
%legend('Test signal','Reconstructed signal','Reconstruction Error');
error = norm(x-xhat,2)^2/norm(x,2)^2;
fprintf('the error occured in the reecovery of signal using FBP is:%d\n',error);