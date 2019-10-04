%%
load('d_08_1_1_4.mat');% Brigde acceleration Data
data = Data(1:20000,2:1:17);
data_mean = mean(data);
[max_mean max_i] = max(data_mean);
[N m] = size(data);
fid = fopen('Input_Output_Details1.txt','w');% Text file storing all the outputs
%%
%%%%%% nosie inclusion %%%%
noise_ind = input('Enter \n 1: if noise present\n 0:if noiseless :');
%%%%%%% compression %%%%%%
method_types = input ('enter \n1: for Douglas Rachford Method \n2: for SPGL1 \n3: for OMP \n4:for CoSaMP\n5: eONEL1 \n6: rONEL1 \n7: for RPCA \n8: Approximate message passing\n9: for gOMP\n10: for Forward Backard Pursuit\n11: for stCSMUSIC\n12:for MUSIC algorithm\n13: for all the above algorithm');
comp_type = input('Enter 1 for Joint Compression and 0 for 0 for Individual compression\n');
%%
if (method_types == 13 || method_types == 7)%%%%  Asking the user Whether to reshape the matrix.
reshape = input('Enter \n 1: if data needs to be reshaped\n 0:if data need not be reshaped\n');
end
if method_types == 13
    method_types = [ 2 3 4 5 6 7 8 9 10 11 12 1];
end
len = length(method_types);

rec_result = cell(len,1);
StdL1 = zeros(len,1);
StdL2 = zeros(len,1);
MinL1 = zeros(len,1);
MaxL1 = zeros(len,1);
MinL2 = zeros(len,1);
MaxL2 = zeros(len,1);
Mean_ErrorL1 = zeros(len,1);
Mean_ErrorL2 = zeros(len,1);
method_names = cell(len,1);
rec_time = zeros(len,1);
figures = cell(len,1);

for i = 1:len    
    method_type = method_types(i);
    fprintf(fid,'\n Method :%d  \n\n',method_type);
tic;
if method_type ~= 7
    [c_noisy_data, K, f_data, s_data,ff_data, A,p] = compressor(data,method_type,comp_type,noise_ind);
else    
    if ~reshape
    [eigv, scores,pp,p_bar,K,pos] = compression_rPCA(data);
    else
    [eigv,scores,pp,p_bar,l] = compression_rPCA_reshape(data);   
    end   
end
c_time = toc;
%%
if method_type ~= 7
    [n m] = size(f_data);
    fprintf(fid,'the size of the actual signal is [ %d X %d ]\n',N,m);
    fprintf(fid,'the size of the comressed signal is [ %d X %d ]\n',p,m);
    fprintf(fid,'the number of peaks in the signal is %d \n',K);
    fprintf(fid,'the time taken for compressionis:%d\n\n',c_time);
else
    if ~reshape
    fprintf(fid,'the size of the actual signal is [ %d X %d ]\n',N,m);
    fprintf(fid,'the size of the comressed signal is [ %d X %d ]\n',K,pp);    
    else
        fprintf(fid,'the size of the actual signal is [ %d X %d ]\n',N,m);
        fprintf(fid,'the size of the reshaped signal is [ %d X %d ]\n',(N/l),(m*l));
        fprintf(fid,'the size of the comressed signal is [ %d X %d ]\n',N/l,pp);        
    end
end
%%
%%%%%%% recovery %%%%%%%
tic;
switch method_type
    case 1
        x = zeros(n,m);
        ift_x = zeros(n,m);
        for is = 1:m   
        [x(:,is) ift_x(:,is)] = test_l1_constraint_n(c_noisy_data(:,is),n,A);       
        end
         method_name = 'DRite';
    case 2
        tau = 0;
        x0  = [];
        sigma =0;
        options= [];
        x = zeros(n,m);
        ift_x = zeros(n,m);
        for is = 1:m 
        [x(:,is),r,g,info] = spgl1(A,c_noisy_data(:,is),tau,sigma,x0,options);
        ift_x(:,is) = ifft(x(:,is),'symmetric');
        end
        method_name = 'SPGl1';
    case 3
        x = zeros(n,m);
        ift_x = zeros(n,m);
        for is = 1:m  
        [x(:,is) ift_x(:,is)] = OMP_n(A,c_noisy_data(:,is),K);
        end
        method_name = 'OMP  ';
    case 4
        x = zeros(n,m);
        ift_x = zeros(n,m);
        for is = 1:m  
        [x(:,is) ift_x(:,is)] = CoSaMP_n(A,c_noisy_data(:,is),K);
        end
        method_name = 'CoSaMP';
        
    case 5
        x = zeros(n,m);
        ift_x = zeros(n,m);
        for is = 1:m  
        [x(:,is) ift_x(:,is)] = eONEL1_n(A,c_noisy_data(:,is),n);
        end
        method_name = 'eONEL1';  
    case 6
        x = zeros(n,m);
        ift_x = zeros(n,m);
        for is = 1:m  
        [x(:,is) ift_x(:,is)] = rONEL1_n(A,c_noisy_data(:,is),n);
        end
        method_name = 'rONEL1';
    case 7         
        if reshape 
        pos = [];
        K =0;
        [ r_data] = recovery_rPCA(scores,eigv,pos,N/l,m*l,p_bar,K,l,reshape);
        else           
        [r_data ] = recovery_rPCA(scores,eigv,pos,N,m,p_bar,K,0,reshape);
         r_data = ifft(r_data);
        end
        method_name = 'rPCA';
    case 8
        T = 1000;
        tol = 1e-15;
        x = zeros(n,m);       
        for is = 1:m  
        x(:,is) = reconstructAmp_n(A, c_noisy_data(:,is), T, tol, s_data(:,is), 0);
        end
        ift_x = ifft(x);
        method_name = 'AMP  ';
    case 9
        S = 20;
        err = 1e-5;
        x = zeros(n,m);        
        for is = 1:m       
        [x(:,is), support, iteration] = gomp(c_noisy_data(:,is), A, K, S, err);
        end
        ift_x = ifft(x);
        method_name = 'GOMP ';
    case 10
        eps = 1e-7;
        fs = 10;
        br = 7;
        TerMode = 'Err';
        display('Running FBP...');
        x = zeros(n,m);        
        for is = 1:m 
        x(:,is) = FBP(c_noisy_data(:,is), A, fs, br, TerMode, eps, K);
        end
        method_name = 'FBP  ';
        ift_x = ifft(x);
    case 11
        [actSet, x] = SeqCSMUSIC(c_noisy_data, A, K);
        ift_x = ifft(x);
        method_name = 'Joint recovery(seqCS_MUSIC)';
    case 12
        
        [actSet, x] = MUSIC(c_noisy_data, A, K);
        ift_x = ifft(x);
        method_name = 'Joint recovery(MUSIC)';
    otherwise
            disp('unknown method.');
end
r_time  = toc;
rec_time(i,1) = r_time;
figures{i,1} = method_name;
%%
fprintf(fid,'Recovery completed\n');
fprintf(fid,'the time taken for the recovery by the "%s method:" is :%d\n',method_name,r_time);
%%%%%%% end of recovery %%%%%
%%
%%%%%%% Full Data Inverse Fourier Transform %%%%%%%   
if method_type ~= 7
ext_x = zeros(length(data),m);
ext_x(1:length(f_data),1:m) = x(1:length(f_data),1:m);
ext_x((length(f_data)+1:length(ff_data)),1:m) = flipud(conj(x(2:length(f_data)-1,1:m)));
ifft_x = ifft(ext_x);
else
    ifft_x = r_data;
end
%%
%%%%%%%  error calcuation  %%%%%%%
errorl1 = zeros(m,1);
errorl2 = zeros(m,1);
for is = 1:m
errorl1(is,1) = norm(data(:,is) - ifft_x(:,is),2)^2/norm(data(:,is),2)^2;
errorl2(is,1) = norm(data(:,is)- ifft_x(:,is),2)^2/norm(data(:,is),2)^2;
end
%r_error = norm(f_data-x,2)^2/norm(f_data)^2;
%fprintf(fid,'the least squares error in the recovery of data by "%s"is: %d\n',method_name,sum(errorl2)/m);
%fprintf(fid,'the least squares error in the recovery of sparse data by "%s"is: %d\n',method_name,r_error/m);
%fprintf(fid,'the l1-norm error in the recovery of data by "%s"is: %d\n\n',method_name,sum(errorl1)/m);

Mean_ErrorL1(i) = mean(errorl1);
Mean_ErrorL2(i) = mean(errorl2);
StdL1(i) = std(errorl1);
StdL2(i) = std(errorl2);
MinL1(i) = min(errorl1);
MaxL1(i) = max(errorl1);
MinL2(i) = min(errorl2);
MaxL2(i)= max(errorl2);
method_names{i,1} = method_name;
xlswrite('ErrorsVsMethods.xlsx',method_names,'A2:A12');
xlswrite('ErrorsVsMethods.xlsx',Mean_ErrorL1,'B2:B12');
xlswrite('ErrorsVsMethods.xlsx',MinL1,'C2:C12');
xlswrite('ErrorsVsMethods.xlsx',MaxL1,'D2:D12');
xlswrite('ErrorsVsMethods.xlsx',StdL1,'E2:E12');
xlswrite('ErrorsVsMethods.xlsx',Mean_ErrorL2,'G2:G12');
xlswrite('ErrorsVsMethods.xlsx',MinL2,'H2:H12');
xlswrite('ErrorsVsMethods.xlsx',MaxL2,'I2:I12');
xlswrite('ErrorsVsMethods.xlsx',StdL2,'J2:J12');
xlswrite('ErrorsVsMethods.xlsx',rec_time,'K2:K12');
figure;
plot(1:1000,abs(ifft_x(1:1000,max_i)));
hgsave(figures{i,1});
rec_result{i,1} = ifft_x ; 

end
save Result.mat rec_result;
fprintf(fid,'Method Name\tL1-Error(avg)\tL2-Error(avg)\n');
for i = 1:len
%fprintf(fid,'Average L1Error occured in method:%d is %d\n',method_type,ErrorL1(i));
%fprintf(fid,'Average L2Error occured in method:%d is %d\n\n',method_type,ErrorL2(i));
fprintf(fid,'%s      \t%.4d     \t%.4d\n',method_names{i},Mean_ErrorL1(i),Mean_ErrorL2(i));
end