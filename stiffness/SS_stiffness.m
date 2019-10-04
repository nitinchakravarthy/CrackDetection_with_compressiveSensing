function [K,MM,C,freq] = SS_stiffness(Length,Number,b,d,an)

nodes=Number+1;
L= Length/Number*ones(1,Number);
Area=b*d;
I=b*d^3/12;

%  Material properties
E=3.45e10*ones(1,Number);
Density = 2400;
    


dof_per_node = 3;
DOF =nodes*dof_per_node;

%END CONDITIONS OF THE BEAM
% 0.FREE
% 1.RESTRAINED

%BC=[UX1 UY1 RZ1 UX2 UY2 RZ2]

BC=[1 1 0 1 1 0];
% an = 0.6 + (0.99-0.6).*rand(Number,1);

%GLOBAL STIFFNESS AND MASS MATRIX INITIALIZATION
K=zeros(DOF);
MM=zeros(DOF);
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
     kuu=an(e).*kuu;
     
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
         for bi=-2:3
             K(a+e*3,bi+e*3)=K(a+e*3,bi+e*3)+kuu(a+3,bi+3);
             MM(a+e*3,bi+e*3)=MM(a+e*3,bi+e*3)+muu(a+3,bi+3);
         end
     end
 end
%APPLYING THE BOUNDARY CONDITION
 BDOF=[1 2 3 DOF-2 DOF-1 DOF];
 for e=6:-1:1
      if BC(1,e)==1
          MM(BDOF(1,e),:)=[];
          MM(:,BDOF(1,e))=[];
          K(BDOF(1,e),:)=[];
          K(:,BDOF(1,e))=[];
      end
 end
 % ==================== Damping matrix ==================================
    [~,B] = eig(K,MM);
    freq  = sqrt(diag(B))/6.28;
    MassM=MM(3:3:size(MM,1),3:3:size(MM,2));
    StiffM=K(3:3:size(K,1),3:3:size(K,2));
    
    T1 = 6.28/freq(1,1);
    T2 = 6.28/freq(2,1);
    B1 = sqrt(B(1,1))/6.28;
    B2 = sqrt(B(2,2))/6.28;
   
    Dp1=.02;Dp2=.02;
    AB = inv([1 B1^2; 1 B2^2]) * [2*B1*Dp1; 2*B2*Dp2];
    alp = AB(1); bet = AB(2);
    C = alp*MM + bet*K;