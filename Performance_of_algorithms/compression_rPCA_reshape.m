function [eigv,scores,pp,p_bar,l] = compression_rPCA_reshape(data)
% Here we compress the data
%data = xlsread('data.xls');
[n m] = size(data); 
l=25;  
w = m*l;v= n/l;        % w and v
C = cell(1,l);   % Creastes a Temporary cell "C"
v_acc_bar=zeros(v,w);

for i = 1:l     % Reshape loop
    C{1,i} = data((i-1)*v+1:i*v,:); 
    v_acc_bar(:,(i-1)*m+1:i*m) = C{1,i};
end
%v_acc_bar_end = v_acc_bar(l*m,:);
%[U_bar_svd, S_bar, V_bar] = svd(v_acc_bar);
%U_bar = U_bar_svd;



K= 60;


 [n_bar,p_bar] = size(data);
pp = 20;
% Now we find the principal components using robpca, there is no need to
% use sparsified data in this analysis as we are just finding principal
% components.
if p_bar < n_bar
    
    x = v_acc_bar;
    %pp = round(0.03*size(x,2)) + 1 ;
 
    fprintf('Compressing the data......\n');
    [lambda, eigv, scores] = robpca(x, pp, @mad) ;
     
     fprintf('Done\n');
    
elseif p_bar > n_bar
    display('The number of columns is greater than the number of rows \n');
     
     x = v_acc_bar' ;
   
     fprintf('Compressing the data......\n');
     [lambda, eigv, scores] = robpca(x, pp, @mad) ;
    
     fprintf('Done.\n');
     
elseif p_bar == n_bar
    display('The number of rows and columns are equal');
   
    x = v_acc_bar ;
   
    fprintf('Compressing the data.....\n');
    [lambda, eigv, scores] = robpca(x, pp, @mad);
    
    fprintf('Done.\n');
    
end