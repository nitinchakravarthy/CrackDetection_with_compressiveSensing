% Here we recover the data
fprintf('Initiating the recovery of data....\n');
tic;
x_rec = scores*eigv' ;
if p < n
   rec_data = x_rec ;
elseif p > n
   rec_data = x_rec' ;
elseif p == n
    rec_data = x_rec ;
end
t_rec = toc;
fprintf('Data recovered\n');
fprintf('The time taken to recover the data is %d \n', t_rec);
% Now we calculate the error associated with data recovery
error = norm(rec_data - data,2)/norm(data,2) ;
fprintf('The error in recovery is %e \n', error);