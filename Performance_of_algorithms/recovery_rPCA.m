function [r_data ] = recovery_rPCA(scores,eigv,pos,n,m,p,K,l,reshape)
% Here we recover the data
fprintf('Initiating the recovery of data....\n');
x_rec = scores*eigv' ;
if p < n
   rec_data = x_rec ;
elseif p > n
   rec_data = x_rec ;
elseif p == n
    rec_data = x_rec ;
end
if ~reshape
r_data = zeros(n,m);
for is = 1:m
r_data(pos(1:K),is) = rec_data(:,is);
end
end

if reshape
    r_data = zeros(n*l,m/l);
   C = cell(1,l);
for i = 1:l      % Reshape loop
    C{1,i} = rec_data(:,(i-1)*m/l+1:i*m/l); 
    r_data((i-1)*n+1:i*n,:) = C{1,i};
end

end

%fprintf('Data recovered\n');
%fprintf('The time taken to recover the data is %d \n', t_rec);
% Now we calculate the error associated with data recovery
%error = norm(rec_data - data,2)/norm(data,2) ;
%fprintf('The error in recovery is %e \n', error);