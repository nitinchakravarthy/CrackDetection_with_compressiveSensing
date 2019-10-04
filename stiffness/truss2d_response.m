 [M,K,C,freq]=truss2d(12);
 x=2:2:size(M,1);
 t = 0.0:.001:.5;
dt = t(2)-t(1);
p1=(1000*sin(2*freq(8)*pi*t)+800*sin(2*freq(10)*pi*t));
pos = x(2:2:end);          % loading locations
L = zeros(size(M,1),length(pos));
L(pos,1:length(pos)) = eye(length(pos));
p = [];
for i = 1:length(pos)
   p = [p; p1]; 
end
P = L * p;
[depl_m, vel_m, acc_m] = NewmarkMethod(M,K,C,P,t);