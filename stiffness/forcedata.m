[M,K,C,freq]=truss2d(12);
dt=1/1000;
t=dt:dt:0.6;
for i=1:length(t)
F1(i)=50*(1-cos(10*pi*i*dt))*sin(30*pi*dt*i);
F2(i)=(65*sin(20*pi*i*dt))+(60*sin(80*pi*i*dt))+(55*sin(160*pi*i*dt));
F3(i)=50*(1-cos(20*pi*i*dt))*sin(40*pi*dt*i);
F4(i)=(50*sin(8*pi*i*dt))+(55*sin(70*pi*i*dt))+(50*sin(150*pi*i*dt));
end
P=zeros(49,length(t));
P(4,:)=F1;
P(16,:)=F2;
P(28,:)=F3;
P(40,:)=F4;
[depl, vel, acc] = NewmarkMethod(M,K,C,P,t);
[f, depl, vel, acc] = force(K,M,C,length(t),dt,[zeros(49,length(t))],acc);
error1=f(4,:)-F1;
SSE1=norm(error1)*100/norm(F1);
error2=f(16,:)-F2;
SSE2=norm(error2)*100/norm(F2);
error3=f(28,:)-F3;
SSE3=norm(error3)*100/norm(F3);
error4=f(40,:)-F4;
SSE4=norm(error4)*100/norm(F4);
finalsse=[SSE1;SSE2;SSE3;SSE4];

