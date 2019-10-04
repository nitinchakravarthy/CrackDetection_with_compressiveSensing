%................................................................
function [stiffness, mass, C, freq,force_data]=mindlin(an)
% MATLAB codes for Finite Element Analysis
% problem19vibrations.m
% Mindlin plate in free vibrations
% antonio ferreira 2008

% clear memory
% clear all;colordef white;clf

% materialsS
N = length(an);
% E  = 10920;     poisson = 0.30;  
% thickness=0.1;
% I=thickness^3/12;
% rho=1;
E1  = 19315.23;    E2 = 19315.23;
poisson = 7.591284e-01;  
thickness=2;
I=thickness^3/12;
rho=1.6e-6;
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

%Mesh generation
L  = 300;    
numberElementsX=sqrt(N);
numberElementsY=numberElementsX;
numberElements=numberElementsX*numberElementsY;


[nodeCoordinates, elementNodes]=...
    rectangularmesh(numberElementsX, numberElementsY,L,L);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes; 

% computation of the system stiffness and mass matrices
[stiffness]=...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,thickness,I,an);

[mass]=...
formMassMatrixMindlinQ4(GDof,numberElements,...
elementNodes,numberNodes,nodeCoordinates,thickness,rho,I);



% % boundary conditions 
[prescribedDof,activeDof,fixedNodeW]=...
EssentialBC('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% G=E/2.6;
% % V : mode shape
% % D : frequency
% % 
% numberOfModes=12;
[V,DD] = eig(stiffness(activeDof,activeDof),...
    mass(activeDof,activeDof)); 
% D = diag(sqrt(DD)*L*sqrt(rho/G));
% [D,ii] = sort(D); ii = ii(1:numberOfModes); 
% VV = V(:,ii);
activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);
% NNN=size(activeDofW);
%     
%     VVV(1:numberNodes,1:12)=0;
%     for i=1:12
%         VVV(activeDofW,i)=VV(1:NNN,i);
%     end
    
NN=numberNodes;N=sqrt(NN);
x=linspace(-L,L,numberElementsX+1);
y=linspace(-L,L,numberElementsY+1);
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


