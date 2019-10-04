function [stiffness, mass]=generateMultipleDamagedFeature_plate(data,boundary_conditions,...
                                            Damages_multiple)
  %%
  % EXPANDING DATA
  thickness = data(3);
  E1        = data(4);
  poisson   = data(6);
  N1        = data(8);   %number of elements along length 
  N2        = data(9);   %number of elements along width
  kapa      = data(10); 
  G12       = data(11);
  GDof      = data(12);
  numberElements = N1*N2;
    E = E1*ones(1,numberElements);
    G12 = G12*ones(1,numberElements);
    I=thickness^3/12;
    load nodeCoordinates.mat nodeCoordinates
    load elementNodes.mat elementNodes
    xx=nodeCoordinates(:,1);
    numberNodes=size(xx,1);
    an=ones(N1*N2,1);
%%
%ALLOWING FOR DAMAGE DATA
for DamageCounter = 1: size(Damages_multiple,1)
E(1,Damages_multiple{DamageCounter,1})=Damages_multiple{DamageCounter,2}*E1;
G12(1,Damages_multiple{DamageCounter,1}) = Damages_multiple{DamageCounter,2}*data(11);
end
E1=E;
E2=E;
C_bending_all = cell(1,numberElements);
C_shear_all = cell(1,numberElements);
% matrix C
% bending part
% C_bending=I*E/(1-poisson^2)*...
%     [1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
for C_counter = 1:numberElements
C_bending_all{1,C_counter} = I*[E1(1,C_counter)/(1-poisson^2) poisson*E2(1,C_counter)/(1-poisson^2) 0;...
    poisson*E1(1,C_counter)/(1-poisson^2) E2(1,C_counter)/(1-poisson^2) 0;...
    0 0 G12(1,C_counter)];

C_shear_all{1,C_counter} = kapa*thickness*G12(1,C_counter)*eye(2);
end
% shear part
% C_shear=kapa*thickness*E/2/(1+poisson)*eye(2);



%%
%FORM STIFFNESS AND MASS MATRICES

[stiffness]=...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear_all,...
    C_bending_all,thickness,I,an);

load mass.mat mass;

end

