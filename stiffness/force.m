function [f, x, xdot, x2dot] = force(K,M,C,pts,dt,f_kn,x2dot_m)

dof = length(M);
% Initialisation
% ===================================================================
x2dot(:,1)=x2dot_m(:,1);    x2dot_c(:,1)=zeros(dof,1);
xdot_c(:,1)=zeros(dof,1);   xdot(:,1)=zeros(dof,1);
x_c(:,1)=zeros(dof,1);  x(:,1)=zeros(dof,1);
f(:,1)=zeros(dof,1);

K_eff = K + 2*C/dt + 4*M/(dt^2);
alpha = 4*M/dt + 2*C;
beta = 2*M;


% computation of load time history
% ===================================================================
for i=2:pts
    
    xdot_c(:,i) = xdot_c(:,i-1) + (x2dot_c(:,i-1) + x2dot_m(:,i))*dt/2;
    x_c(:,i) = x_c(:,i-1)+ (xdot_c(:,i-1) + xdot_c(:,i))*dt/2;
    
    f_un(:,i) = M*x2dot_m(:,i) + C*xdot_c(:,i) + K*x_c(:,i)-f_kn(:,i);
    f(:,i) = f_un(:,i)+f_kn(:,i);
    df(:,i-1)=f(:,i)-f(:,i-1)+alpha*xdot_c(:,i-1) + beta*x2dot_m(:,i-1);
    dx(:,i-1) = K_eff\df(:,i-1);
    dxdot(:,i-1) = 2*dx(:,i-1)/dt - 2*xdot_c(:,i-1);
    dx2dot(:,i-1) =  4*dx(:,i-1)/dt^2 -4*xdot_c(:,i-1)/dt -2*x2dot(:,i-1);
    xdot(:,i)= xdot_c(:,i-1) + dxdot(:,i-1);
    x(:,i)= x_c(:,i-1) + dx(:,i-1);
    x2dot(:,i) = x2dot_c(:,i-1) + dx2dot(:,i-1);
    x2dot_c(:,i)=x2dot(:,i);
    xdot_c(:,i) = xdot(:,i);
    x_c(:,i) = x(:,i);
%     ==============================================================
end
