function loading = create_loading_onGrid_load( fmax, l_w, unl_w,pos,mass,N1,N2 )

% Time step
Fs = 10000;                     % Sampling frequency
dt = 1/Fs;                      % Sample time
L = 1000;                       % Length of signal
t = (0:L-1)*dt;                 % Time vector
%{
omega=freq(5:5+numberModes);
  amp=0;
  for count=1:numberModes
      amp=amp+1E1*rand(1,1)*sin(omega(count)*t)...
             +1E1*rand(1,1)*cos(omega(count)*t);
  end
%}
Force = zeros(length(fmax),L);
for i = 1:length(fmax)
Force(i,:) = impact_load(fmax(i), l_w(i), unl_w(i), t, Fs);
end
%{
L = zeros(size(mass,1),length(pos));
  L(pos,1:length(pos)) = eye(length(pos));
  loading = L * repmat(Force,[length(pos),1]);
%}
loading = zeros(size(mass,1),length(Force));
for i = 1:length(fmax)
loading(pos(i),:) = Force(i,:);
end
  %loading = repmat(amp',[1,length(x)]);
  
end

