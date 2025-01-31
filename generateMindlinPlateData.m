function [stiffness, mass, C, freq,force_data]=generateMindlinPlateData(length,width,seed_length)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%Plate Data

N1=length/seed_length; %along length %aligned with Y-axis
N2=width/seed_length; %along width %aligned with X-axis
an=ones(N1*N2,1);
E1  = 19315.23;    
E2 = 19315.23;
poisson = 7.591284e-01;  
thickness=0.25;
I=thickness^3/12;
rho=1.6e-6;
boundary_conditions = 'ssss';
save boundary_conditions.mat boundary_conditions
severityScale=0.5:0.1:0.9;
save severityScale.mat severityScale
% kapa=0.8601; % cccc / cccf case
% kapa=0.822; % scsc case
kapa=5/6;  % ssss case 
G12 = 36832.12;

% matrix C
% bending part
% C_bending=I*E/(1-poisson^2)*...
%     [1 poisson 0;poisson 1 0;0 0 (1-poisson)/2]; 
C_bending=I*[E1/(1-poisson^2) poisson*E2/(1-poisson^2) 0;...
    poisson*E1/(1-poisson^2) E2/(1-poisson^2) 0;...
    0 0 G12];

% shear part
% C_shear=kapa*thickness*E/2/(1+poisson)*eye(2);
C_shear=kapa*thickness*G12*eye(2);

% Saving Data
data   =zeros(20,1);
  data(1)=length;
  data(2)=width;
  data(3)=thickness;
  data(4)=E1;
  data(5)=E2;
  data(6)=poisson;
  data(7)=rho;
  data(8)=N1;   %number of elements along length 
  data(9)=N2;   %number of elements along width
  data(10)=kapa; 
  data(11)= G12;
  
 %Mesh generation
L  = length;   %length
W = width;
numberElementsX=N2;
numberElementsY=N1;
numberElements=numberElementsX*numberElementsY;


[nodeCoordinates, elementNodes]=...
    rectangularmesh(numberElementsX, numberElementsY,W,L);
save nodeCoordinates.mat nodeCoordinates
save elementNodes.mat elementNodes
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes; 
data(12) = GDof;
save data.mat data
% computation of the system stiffness and mass matrices
[stiffness]=...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,thickness,I,an);

save stifness.mat stiffness
[mass]=...
formMassMatrixMindlinQ4(GDof,numberElements,...
elementNodes,numberNodes,nodeCoordinates,thickness,rho,I);

save mass.mat mass

% % boundary conditions 
[prescribedDof,activeDof,fixedNodeW]=...
EssentialBC(boundary_conditions,GDof,xx,yy,nodeCoordinates,numberNodes);

save prescribedDof.mat prescribedDof
save activeDof.mat activeDof
save fixedNodeW.mat fixedNodeW

[V,DD] = eig(stiffness(activeDof,activeDof),...
    mass(activeDof,activeDof)); 
save V.mat V
save DD.mat DD

activeDofW=setdiff(1:numberNodes',fixedNodeW);
stiffness=stiffness(activeDof,activeDof);
mass=mass(activeDof,activeDof);

xxx=nodeCoordinates(activeDofW,1);
yyy=nodeCoordinates(activeDofW,2);
% drawing Eigenmodes
% drawEigenmodes2D(x,y,VVV,NN,N,D)

B1 = sqrt(DD(1,1))/6.28;
B2 = sqrt(DD(2,2))/6.28;

Dp1=.02;Dp2=.02;
AB = inv([1 B1^2; 1 B2^2]) * [2*B1*Dp1; 2*B2*Dp2];
alp = AB(1); bet = AB(2);
C = alp*mass + bet*stiffness;
freq = sqrt (DD)/6.28;
freq=(sort(diag(freq)));
force_data = struct('d1',{GDof},'d2',{numberElements},'d3',...
    {elementNodes},'d4',{numberNodes},'d5',{nodeCoordinates},...
    'd6',{activeDof});

