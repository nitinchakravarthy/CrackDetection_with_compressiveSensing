function [z,method_name]=evaluateOptimumZ(A,vDamage,method_types)
% A              =  corresponding Dictionary
% vDamage_group =  Signal responses of a particular sensor for different 
%                   damage locations and magnitudes 
% method_types   =  Which method to be used as the case numbers
%{
 B = A;
A = [];
v = vDamage_group;
vDamage_group = [];
for i = 1:2:size(B,1)-1
    A = [ A;    B(i,:) ];
    vDamage_group = [vDamage_group; v(i,:)];
end
%}

 
len = length(method_types);
for i = 1:len    
    method_type = method_types(i);
   
switch method_type   
    case 1    
      zo = A'*vDamage;
        z = l1qc_logbarrier(zo, A, [], vDamage, 0.001);
        
        method_name ='l1-magic';
    
    case 2  
        n=size(A,2);
        z = test_l1_constraint_n(vDamage,n,A);       
        
        method_name = 'DRite';
    case 3
        tau = 0;
        sigma =0;
        options= [];
        zo = A'*vDamage;
        z = spgl1(A,vDamage,tau,sigma,zo,options); 
        
        method_name = 'SPGl1';
    case 4
        K=10; 
        z = OMP( A,vDamage, K,[], [] );
        %z= OMP_n(A,vDamage,K);
        
        method_name = 'OMP  ';
    case 5
        K =10;
      z = CoSaMP_n(A,vDamage,K);
       
        method_name = 'CoSaMP';

    case 6
        S = 1;
        err = 1e-5;
        K=20;
        z = gomp(vDamage, A, K, S, err);
     
        method_name = 'GOMP ';
    case 7
        eps = 1e-7;
        fs = 10;
        br = 7;
        TerMode = 'Err';
        K=20;
        z = FBP(vDamage, A, fs, br, TerMode, eps, K);
        
        method_name = 'FBP  ';
    
    case 8
        K = 20;
        z = SeqCSMUSIC(vDamage, A, K);
        
        method_name = 'Joint recovery(seqCS_MUSIC)';
    case 9
        K = 20;
        z = SOMP(vDamage, A, K);
        
        method_name = 'Joint recovery(MUSIC)';
    case 10 
        z = l1_homotopy_rni_load(real(vDamage),real(A));
        method_name = 'L1 homotopy';
    case 11
        K = 20;
         [z,~] = NNS_greedy(real(A), real(vDamage), K,'alg','sp','nonneg','on', 'simult','off','maxIter',1000);
         method_name = 'NNS_SP';
    case 12
        K = 20;
        [z,~] = NNS_greedy(real(A), real(vDamage), K,'alg','omp','nonneg','on', 'simult','off','maxIter',1000);
         method_name = 'NNS_OMP';
    case 13
        K = 20;
        [z,~] = NNS_greedy(real(A), real(vDamage), K,'alg','htp','nonneg','on', 'simult','off','maxIter',1000);
         method_name = 'NNS_HTP';
    case 14
        K = 20;
        [z,~] = NNS_greedy(real(A), real(vDamage), K,'alg','cosamp','nonneg','on', 'simult','off','maxIter',1000);
          method_name = 'NNS_COSAMP';
    otherwise
            disp('unknown method.');
end
end
end