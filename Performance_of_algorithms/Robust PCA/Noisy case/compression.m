% Here we compress the data
%data = xlsread('data.xls');
data = c_data;
[n,p] = size(data);
noise = randn(n,p);
noisy = data + noise ;
t = (1 : 1: n);
% Now we find the principal components using robpca, there is no need to
% use sparsified data in this analysis as we are just finding principal
% components.
if p < n
    tic;
    x = noisy;
    %pp = round(0.03*size(x,2)) + 1 ;
    pp = 5;
    fprintf('Compressing the data......\n');
    [lambda, eigv, scores] = robpca(x, pp, @mad) ;
     t_comp = toc ;
     fprintf('Done\n');
     fprintf('The time taken for the compression is %d \n', t_comp);
elseif p > n
    display('The number of columns is greater than the number of rows \n');
     tic;
     x = noisy' ;
     pp = round(0.03*size(x,2)) ;
     fprintf('Compressing the data......\n');
     [lambda, eigv, scores] = robpca(x, pp, @mad) ;
     t_comp = toc ;
     fprintf('Done.\n');
     fprintf('The time taken for compression is %d \n',t_comp);
elseif p == n
    display('The number of rows and columns are equal');
    tic;
    x = noisy ;
    pp = round(0.03*size(x,2));
    fprintf('Compressing the data.....\n');
    [lambda, eigv, scores] = robpca(x, pp, @mad);
    t_comp = toc;
    fprintf('Done.\n');
    fprintf('The time taken for the compression is %d \n', t_comp);
end