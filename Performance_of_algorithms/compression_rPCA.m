function [eigv, scores,pp,p_bar,K,pos] = compression_rPCA(data)
% Here we compress the data
%data = xlsread('data.xls');
K= 300;
[n m] = size(data); 
f_data = fft(data);

s_data = zeros(K,m);

[~,pos] = sort(f_data(:,1),'descend');
s_data(:,1) = f_data(pos(1:K),1);
 for is = 2:m
 s_data(:,is) = f_data(pos(1:K),is);
 end
[n_bar,p_bar] = size(s_data);
t = (1 : 1: n_bar);
% Now we find the principal components using robpca, there is no need to
% use sparsified data in this analysis as we are just finding principal
% components.
if p_bar < n_bar
    tic;
    x = s_data;
    %pp = round(0.03*size(x,2)) + 1 ;
    pp = 16;
    fprintf('Compressing the data......\n');
    [lambda, eigv, scores] = robpca(x, pp, @mad) ;
     t_comp = toc ;
     fprintf('Done\n');
     fprintf('The time taken for the compression is %d \n', t_comp);
elseif p_bar > n_bar
    display('The number of columns is greater than the number of rows \n');
     tic;
     x = s_data' ;
     pp = round(0.03*size(x,2)) ;
     fprintf('Compressing the data......\n');
     [lambda, eigv, scores] = robpca(x, pp, @mad) ;
     t_comp = toc ;
     fprintf('Done.\n');
     fprintf('The time taken for compression is %d \n',t_comp);
elseif p_bar == n_bar
    display('The number of rows and columns are equal');
    tic;
    x = s_data ;
    pp = round(0.03*size(x,2));
    fprintf('Compressing the data.....\n');
    [lambda, eigv, scores] = robpca(x, pp, @mad);
    t_comp = toc;
    fprintf('Done.\n');
    fprintf('The time taken for the compression is %d \n', t_comp);
end