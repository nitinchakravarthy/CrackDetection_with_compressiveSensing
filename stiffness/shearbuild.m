
function [M,K,C,freq]=shearbuild(an)

storey = 20;
% storey_damaged = input('Enter the damaged storey : ');      %   [11]
% storey_damaged=5;
% Temp_Vari = input('Enter the temperature : ');
% stiff_reduc = input('Enter the % of stiffness reduction : ');
stiff_reduc =0;
mat = 25;
fact=concrete(mat);
Factor = fact;

% ========================= Mass matrix =============================
M = zeros(storey,storey);
for i = 1:10
    M(i,i) = 2000;   % mass of first 10 stories in kgm
end
for i = 11:storey
    M(i,i) = 1500;   % mass of last 10 stories in kgm
end

%========================= Stiffness matrix =============================
%Find k for individual members
kk = zeros(storey,1);

for i = 1:10
    kk(i,1) = Factor*50000e3*an(i); % stiffness of 1-10 stories in N/m
end
for i = 11:15
    kk(i,1) = Factor*40000e3*an(i); % stiffness of 11-15 stories in N/m
end 
for i = 16:storey      
    kk(i,1) = Factor*30000e3*an(i); % stiffness of 16-20 stories in N/m
end

% for n = storey_damaged
%     kk(n,1) = ((100-stiff_reduc)/100*kk(n,1));
% end

%======== Assembling Of Stiffness Matrix ========== %
K = zeros(storey,storey);
K(1,1) = kk(1)+kk(2); 
K(storey,storey) = kk(storey);
K(storey,(storey-1)) = -kk(storey);
K((storey-1),storey) = -kk(storey);
for a = (storey-1):-1:2
    K(a,a) = kk(a)+kk(a+1);
    K(a,a-1) = -kk(a);
    K(a-1,a) = -kk(a);
end
%================== Daming matrix ================= %

[~, D]=eig(K,M);
omg=diag(sqrt(D));
tmp=[omg(1)^2 1;omg(2)^2 1];tmp1 = [2*omg(1)*.01; 2*omg(2)*0.01];
tmp3 = tmp\tmp1;
C=tmp3(1)*K +tmp3(2)*M; %raleigh damping
clear tmp tmp1 tmp3;
freq=omg./(2*pi);


