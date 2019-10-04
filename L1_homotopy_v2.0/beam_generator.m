%%%%%%%%%%%% Random Force Generator Code) %%%%%%%%%%%
%rng(0);
freq1= 100;
freq2=1000;
fs=1000;
amp=1;
lenS=1;
t=1/fs:1/fs:lenS; %Generates a time discretization vector
disp(length(t));
f= randperm(1000)';
f = f(1:10,1);
wwf = f;
disp(length(f));
force =zeros(1000,1);
for i=1:length(f)
force = force + amp*(cos((2*pi*f(i,1)*t')) );
end
save force.mat force


%%%%%%%%%%%%%%%%%%%%%%%%
% Beam stiffness, mass and damping matrices %
%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('InputDetails.txt','a');
%Y=input('Enter the Material Of the Beam 1-concrete 2- Steel ===>      ');
Y=2;
%Length= input('Enter the Length Of the Beam(m)===>    ');
Length=1;
%Number= input('Enter the number Of Elements===>       ');
Number=10;
nodes=Number+1;
L= Length/Number*ones(1,Number);  
%b=input('Enter the Width Of the Beam(m)===>      ');
%d=input('Enter the Depth Of the Beam(m)===>      ');
%pts=input('enter number of time steps : ');
%dt=input('enter the length of a time step : ');
b=0.014;
d=0.014;
pts=5000;
dt=1/10000;
time=dt:dt:pts*dt;
Area=b*d;
I=b*d^3/12;
fprintf(fid,'Type of Beam : Cantilever\n');
if Y==1
    E=2.5e10*ones(1,Number);
    Density = 2400;
    fprintf(fid,'Concrete Material\n');
end
if Y==2
    E = 2.0e11*ones(1,Number);
    Density = 7800;
    fprintf(fid,'Steel Material\n');
end
fprintf(fid,'Beam dimensions : %4.3f * %4.3f\n',b,d);
fprintf(fid,'Beam length : %5.3f\n',Length);
fprintf(fid,'Time Steps : %5d\n',pts);
fprintf(fid,'Length of a time step : %5.4f\n',dt);
dof_per_node = 3;
DOF =nodes*dof_per_node;

%END CONDITIONS OF THE BEAM     
% 0.FREE
% 1.RESTRAINED 

%BC=[UX1 UY1 RZ1 UX2 UY2 RZ2]

 BC=[1 1 1 0 0 0];                                                                                      

%GLOBAL STIFFNESS AND MASS MATRIX INITIALIZATION
 K=zeros(DOF);
 M=zeros(DOF);
%ASSEMBLING THE MATRICES
 for e=1:Number
%ELEMENT STIFFNESS MATRIX
     kuu=zeros(6);
     kuu(1,1)=E(e)*Area/L(e);
     kuu(1,4)=-E(e)*Area/L(e);
     kuu(2,2)=(12*E(e)*I)/(L(e)^3);
     kuu(2,3)=(6*E(e)*I)/(L(e)^2);
     kuu(2,5)=-(12*E(e)*I)/(L(e)^3);
     kuu(2,6)=(6*E(e)*I)/(L(e)^2);
     kuu(3,2)=(6*E(e)*I)/(L(e)^2);
     kuu(3,3)=4*E(e)*I/L(e);
     kuu(3,5)=-(6*E(e)*I)/(L(e)^2);
     kuu(3,6)=2*E(e)*I/L(e);
     kuu(4,1)=-E(e)*Area/L(e);
     kuu(4,4)=E(e)*Area/L(e);
     kuu(5,2)=-(12*E(e)*I)/(L(e)^3);
     kuu(5,3)=-(6*E(e)*I)/(L(e)^2);
     kuu(5,5)=(12*E(e)*I)/(L(e)^3);
     kuu(5,6)=-(6*E(e)*I)/(L(e)^2);
     kuu(6,2)=(6*E(e)*I)/(L(e)^2);
     kuu(6,3)=2*E(e)*I/L(e);
     kuu(6,5)=-(6*E(e)*I)/(L(e)^2);
     kuu(6,6)=4*E(e)*I/L(e);
%ELEMENT MASS MATRIX
     muu=zeros(6);
     muu(1,1)=140;
     muu(1,4)=70;
     muu(2,2)=156;
     muu(2,3)=22*L(e);
     muu(2,5)=54;
     muu(2,6)=-13*L(e);
     muu(3,2)=22*L(e);
     muu(3,3)=4*L(e)^2;
     muu(3,5)=13*L(e);
     muu(3,6)=-3*L(e)^2;
     muu(4,1)=70;
     muu(4,4)=140;
     muu(5,2)=54;
     muu(5,3)=13*L(e);
     muu(5,5)=156;
     muu(5,6)=-22*L(e);
     muu(6,2)=-13*L(e);
     muu(6,3)=-3*L(e)^2;
     muu(6,5)=-22*L(e);
     muu(6,6)=4*L(e)^2;
     muu=muu*Density*Area*L(e)/420;
     for a=-2:3
         for b=-2:3
             K(a+e*3,b+e*3)=K(a+e*3,b+e*3)+kuu(a+3,b+3);
             M(a+e*3,b+e*3)=M(a+e*3,b+e*3)+muu(a+3,b+3);
         end
     end
 end
%APPLYING THE BOUNDARY CONDITION
 BDOF=[1 2 3 DOF-2 DOF-1 DOF];
 for e=6:-1:1
      if BC(1,e)==1
          M(BDOF(1,e),:)=[];
          M(:,BDOF(1,e))=[];
          K(BDOF(1,e),:)=[];
          K(:,BDOF(1,e))=[];
      end
 end
save massmatrix.mat M;
save stiffnessmatrix.mat K;
% ==================== Damping matrix ==================================
[eigen_v_disp,B] = eig(K,M);
freq  = sqrt(diag(B))/6.28;
Fundamental_Frequency=freq(1,1);
fprintf(fid,'Fundamental Frequency : %5.4f\n',Fundamental_Frequency);
Second_Frequency=freq(2,1);
fprintf(fid,'Second Frequency : %2.4f\n',Second_Frequency);
ffreq = min(freq)*0.5;
%%xlswrite('modeshapes17_2Cant.xlsx',A);
%%xlswrite('frequencies17_2Cant.xlsx',freq);
MassM=M(2:3:size(M,1),2:3:size(M,2));
StiffM=K(2:3:size(K,1),2:3:size(K,2));
FE_modeshape=eigen_v_disp(2:3:size(eigen_v_disp,1),1:ceil(size(eigen_v_disp,2)/2));
FE_freqs=freq;

T1 = 6.28/freq(1,1);
T2 = 6.28/freq(2,1);
B1 = sqrt(B(1,1))/6.28;
B2 = sqrt(B(2,2))/6.28;

Dp(1) =0.01;
Dp(2) =0.005;
fprintf(fid,'Damping values : %2.4f & %2.4f\n\n\n',Dp(1),Dp(2));
AB = inv([1 B1^2; 1 B2^2]) * [2*B1*Dp(1); 2*B2*Dp(2)];
alp = AB(1); bet = AB(2);
disp(['Rayleigh coefficients: alpha =  ' num2str(alp) '   beta = ' num2str(bet)]);
C = alp*M + bet*K;
save dampingmatrix.mat C;
load('stiffnessmatrix.mat');
load('dampingmatrix.mat');
load('massmatrix.mat');


%%%%%%%%%%%%%%%%%%%%%
%Newmark generator for time history responses %
%%%%%%%%%%%%%%%%%%%%%
cidep = zeros(30,1);
civit = zeros(30,1);
tf=1;
pas = 1e-3;
pf = [zeros(28,1000);force';zeros(1,1000)];
taille=length([pas:pas:tf]');
nddl=length(M(:,1));
gamma=0.5;beta=0.25;eps=10^-7;
%alpha=1e-8;gamma=0.5+alpha;beta=0.25*(1+alpha)^2;eps=10^-7;
if isempty(C)==1
    C=zeros(nddl,nddl);
end
    
    dep=zeros(nddl,taille);vit=zeros(nddl,taille);acc=zeros(nddl,taille);
    
    S=M+gamma*pas*C+beta*pas*pas*K;invS=inv(S);
    A1=(1-gamma)*pas;A2=(0.5-beta)*pas*pas;A3=pas*gamma;A4=pas*pas*beta;A5=taille;
    
    acc(:,1)=M\(pf(:,1)-C*civit-K*cidep);
    vit(:,1)=civit;
    dep(:,1)=cidep;
    
    for iter=2:A5
        
        % Prediction
        
        vit(:,iter)=vit(:,iter-1)+A1*(acc(:,iter-1));
        dep(:,iter)=dep(:,iter-1)+pas*vit(:,iter-1)+A2*acc(:,iter-1);
        
        % Computation of the acceleration
        
        acc(:,iter)=invS*(pf(:,iter)-C*vit(:,iter)-K*dep(:,iter));
        
        % Correction
        
        vit(:,iter)=vit(:,iter)+A3*acc(:,iter);
        dep(:,iter)=dep(:,iter)+A4*acc(:,iter);
    end
    vertcl= 2:3:29;
%%this is for taking only vertical degrees of freedom responses
%%applicable for cantilever beam not for shear building and    %%spring mass system
v_disp = zeros(10,iter);
v_vel = zeros(10,iter);
v_acc = zeros(10,iter);
for i= 1:10
    v_disp(i,1:iter)= dep(vertcl(i),1:iter); 
    v_vel(i,1:iter)= vit(vertcl(i),1:iter);
    v_acc(i,1:iter)= acc(vertcl(i),1:iter); 
end
v_disp_end = v_disp(10,1:iter);
v_vel_end = v_vel(10,1:iter);
v_acc_end = v_acc(10,1:iter);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
