% NEWMARK'S METHOD : LINEAR SYSTEM
% Variable Description :
% INPUT :
%       M - Mass Matrix (in modal coordinates)
%       K - Stiffness Matrix (in modal coordinates)
%       C - Damping Matrix (in modal coordinates)
%       P - Force Matrix (in modal coordinates)
%       acceleration - Type of Newmark's Method to be used 
%       
% OUTPUT :
%        depl - modal displacement's 
%        vel - modal velocities
%        accl - modal accelerations
%        t - time values at which integration is done
%--------------------------------------------------------------------------
function [depl, vel, acc] = NewmarkMethod(M,K,C,P,t)

gaama = 1/2 ;beta = 1/4 ; % Newmarks average acceleration

% Time step
dt = t(2)-t(1);
nt=length(t);
n = length(M);

% Constants used in Newmark's integration
a1 = gaama/(beta*dt) ;      a2 = 1/(beta*dt^2) ;
a3 = 1/(beta*dt) ;          a4 = gaama/beta ;
a5 = 1/(2*beta) ;           a6 = (gaama/(2*beta)-1)*dt ;


depl = zeros(n,nt) ; % m
vel = zeros(n,nt) ; % m/sec
acc = zeros(n,nt) ; % m/sec^2
% Initial Conditions
depl(:,1) = zeros ; % m
vel(:,1) = zeros ; % m/sec
acc(:,1) = M\(( P(:,1))-(C*vel(:,1))-(K*depl(:,1))) ;

Kcap = (K)+(a1*C)+a2*M ;

a = a3*M+a4*C ;
b = a5*M+a6*C ;

% Time step starts
for i = 1:nt-1
    delP = (P(:,i+1)-P(:,i)) +a*vel(:,i)+b*acc(:,i) ;
    delu = Kcap\delP ;
    delv = a1*delu-a4*vel(:,i)-a6*acc(:,i) ;
    dela = a2*delu-a3*vel(:,i)-a5*acc(:,i);
    depl(:,i+1) = depl(:,i)+delu ;
    vel(:,i+1) = vel(:,i)+delv ;
    acc(:,i+1) = acc(:,i)+dela ; 
end

