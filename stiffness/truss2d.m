% MATLAB code for Finite Element Analysis of truss
function [M,K,C,freq]=truss2d(n)
% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E=2.06*10^11; A=0.0016; EA=E*A; rho=A*7850;

% generation of coordinates and connectivities
for i=1:n+1
    nodecoo(i*2-1,1)=0+(i-1)*1;
    nodecoo(i*2,1)=0+(i-1)*1;
    nodecoo(i*2-1,2)=0;
    nodecoo(i*2,2)=1;
end
for i=1:n
    elnodes(i*5-4,1)=1+(i-1)*2;elnodes(i*5-4,2)=i*2;
    elnodes(i*5-3,1)=1+(i-1)*2;elnodes(i*5-3,2)=i*2+1;
    elnodes(i*5-2,1)=1+(i-1)*2;elnodes(i*5-2,2)=i*2+2;
    elnodes(i*5-1,1)=(i)*2; elnodes(i*5-1,2)=i*2+1;
    elnodes(i*5,1)=(i)*2;elnodes(i*5,2)=i*2+2;
end
elnodes(n*5+1,1)=n*2+1;elnodes(n*5+1,2)=n*2+2;
numel=size(elnodes,1);
numnode=size(nodecoo,1);
xx=nodecoo(:,1);
yy=nodecoo(:,2);
gdof=2*numnode;
U=zeros(gdof,1);
force=zeros(gdof,1);

% computation of the system stiffness & Mass matrix
K=zeros(gdof); 
M=zeros(gdof); 
% computation of the system stiffness matrix
for e=1:numel; 
  % eldof: element degrees of freedom (Dof)
  index=elnodes(e,:)   ;       
  eldof=[ index(1)*2-1 index(1)*2 index(2)*2-1 index(2)*2] ;
  xa=xx(index(2))-xx(index(1));
  ya=yy(index(2))-yy(index(1));
  l=sqrt(xa*xa+ya*ya); % length of the element
  C=xa/l;
  S=ya/l;   
    k1=EA/l*[C*C C*S -C*C -C*S; C*S S*S -C*S -S*S;
        -C*C -C*S C*C C*S;-C*S -S*S C*S S*S];    
  K(eldof,eldof)= K(eldof,eldof)+k1;
  m1=rho*l*100/2*[1 0 0 0; 0 1 0 0; 0 0 1 0 ; 0 0 0 1;];
  M(eldof,eldof)= M(eldof,eldof)+m1;
end 
% boundary conditions and solution
prescribedDof=[1 2 (gdof-2)]';%(n+1)*4-2
activedof=setdiff([1:gdof]',[prescribedDof]);
K=K(activedof,activedof);
M=M(activedof,activedof);
% drawingMesh(nodecoo,elnodes,'L2','k');
% ===============================================================
[~, D]=eig(K,M);
omg=diag(sqrt(D));
tmp=[omg(1)^2 1;omg(2)^2 1];tmp1 = [2*omg(1)*.0055; 2*omg(2)*0.0039];
tmp3 = tmp\tmp1;
C=tmp3(1)*K +tmp3(2)*M; %raleigh damping
% C=.1523*M+(4.6503*10^-5)*K;
clear tmp tmp1 tmp3;
omega=diag(sqrt(D));
freq=omega/6.28;
% ===============================================================


