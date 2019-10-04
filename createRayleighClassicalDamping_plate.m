function [stiffness,mass,damping]=createRayleighClassicalDamping_plate(stiffness,mass,activeDof)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION
% Creates Rayleigh Classical Damping (partitioned to active DOFs) from full-stiffness and full-damping
% matrices
%
% INPUTS
% stiffness              : full, contains all DOFS, even the fixed ones
% mass                   : full, contains all DOFS, even the fixed ones
% activeDof              : a vector containing numbers of not fixed DOFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Damping matrix needs to be created only for the partitioned K and M
% matrix
stiffness=stiffness(activeDof,activeDof);
mass=mass(activeDof,activeDof);
  [V,D] = eig(stiffness,mass); 
  D = diag(D);
  B1 = sqrt(D(1))/6.28;
  B2 = sqrt(D(2))/6.28;
  Dp(1)=0.01;                                                                % modal damping for first  state
  Dp(2)=0.01;                                                                % modal damping for second state
  AB = inv([1 B1^2; 1 B2^2]) * [2*B1*Dp(1); 2*B2*Dp(2)];
  alp = AB(1); bet = AB(2);
  damping = alp*mass + bet*stiffness;                                                       % already partitionedend
end  