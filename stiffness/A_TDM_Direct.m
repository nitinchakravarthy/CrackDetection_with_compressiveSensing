% Clear all existing data
clear all;clc;close all;

%% Input as FEM models
% ===================================================================
% SIMPLY SUPPORTED BEAM MODEL
% Function format: [K,M,C,freq] = SS_stiffness(Length,N,b,d,theta0);
% N=20;                % Number of elements
% theta0=ones(N,1); % design variables
% [K,M,C,freq] = SS_stiffness(3,N,0.3,0.5,theta0);
% x=3:3:size(M,1);
% -------------------------------------------------------------------
% SHEAR BUILDING
% [M,K,C,freq]=shearbuild(ones(20,1));  % 20 storey shear building of concrete 
% x=1:size(M,1); 
%     
% -------------------------------------------------------------------
% TRUSS BRIDGE
%function call format: [M,K,C]=truss2d(an); an=number of bays.
 [M,K,C,freq]=truss2d(12);
 x=2:2:size(M,1);
% -------------------------------------------------------------------
% PLATE 
% function call format: [M,K,C]=truss2d(an); an=number of bays.
%N1=10;N2=10;
%[K,M,C,freq,force_data]=mindlin(ones(N1*N2,1));
%x = 1:((N1-1)*(N2-1));
% ===================================================================

% Time step
t = 0.0:.001:.5;
dt = t(2)-t(1);
%% Loading
% -------------------------------------------------------------------
% r=rand(150,1);w=randi([10 160],150,1);p=zeros(1,length(t));
% for i=1:length(r);
%     p=p+sin(w(i)*t)*r(i)*10^5;
% end
p1=(1000*sin(2*freq(8)*pi*t)+800*sin(2*freq(10)*pi*t));
p2=(-1000*sin(2*freq(15)*pi*t)+800*sin(2*freq(5)*pi*t));

pos = x([30 60]);          % loading locations
L = zeros(size(M,1),length(pos));
L(pos,1:length(pos)) = eye(length(pos));
P = L * [p1; p2;];
%% To obtain measurement data for load identification, Newmarks-beta
% algorithm is used
[depl_m, vel_m, acc_m] = NewmarkMethod(M,K,C,P,t);

%% Guyans Reduction
y=setdiff(1:size(M,1),x);
acc_m=acc_m(x,:); % measured vertical acceleration data
[M, K, C, f_kn] = guyan(M, K, C, P, length(t),x,y);
dof = length(M);
f_kn = zeros(size(f_kn));
% 
%% Force Identification

[f, depl, vel, acc] = force(K,M,C,length(t),dt,f_kn,acc_m);

%% Error plots

% SSE1=norm(abs(P(DOF(1),:)-f(DOF(1),:)))*100/norm(abs(P(DOF(1),:)));
% SSE2=norm(abs(P(DOF(2),:)-f(DOF(2),:)))*100/norm(abs(P(DOF(2),:)));
% 
% 
% plot(t,P(DOF(1),:)/10^7,'r-','linewidth',3);hold all;
% plot(t(1:5:501),f(DOF(1),(1:5:501))/10^7,'b-o','linewidth',1); hold all;
% grid('on');
% % plot(t,p_bar(:,1),'b','linewidth',1); hold all;
% xlabel('Time (sec)','Fontsize',12); ylabel('Load (kN)','Fontsize',12);
% legend('Actual loading','Identified');
% axis([0 .5 -1500 2000]);set(gca,'FontSize', 12);
% figure(2);
% plot(t,P(DOF(2),:)/10^7,'r-','linewidth',3);hold all;
% plot(t(1:5:501),f(DOF(2),(1:5:501))/10^7,'b-o','linewidth',1); hold all;
% grid('on');
% % plot(t,p_bar(:,2),'b','linewidth',1); hold all;
% xlabel('Time (sec)','Fontsize',12); ylabel('Load (kN)','Fontsize',12);
% legend('Actual loading','Identified');
% axis([0 .5 -1500 2000]);set(gca,'FontSize', 12);





