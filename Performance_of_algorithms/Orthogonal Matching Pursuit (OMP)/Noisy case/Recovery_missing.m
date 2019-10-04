data =v_acc';
[n m] = size(data); 
l=100;
w = m*l;v= n/l;        % w and v
C = cell(1,l);   % Creastes a Temporary cell "C"
v_acc_bar=zeros(v,w);

for i = 1:l     % Reshape loop
    C{1,i} = data((i-1)*v+1:i*v,:); 
    v_acc_bar(:,(i-1)*m+1:i*m) = C{1,i};
end
A = v_acc_bar;
B = rand(size(A))<.9;
lamnbda_tol = 10;
tol = 1e-7;
N = 100;
fprintf('Completion using nuclear norm minimization... \n');
[CompletedMat, ier] = MatrixCompletion(A.*B, B,N, 'nuclear', lamnbda_tol, tol, 0);




