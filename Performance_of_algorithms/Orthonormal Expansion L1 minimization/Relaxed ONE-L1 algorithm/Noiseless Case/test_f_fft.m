i = 1;
new_data = fft(force);
K  = 20;
[~, pos] = sort(new_data(:,i),'descend');
    s_data(pos(1:K),i) = new_data(pos(1:K),i) ;
    %figure;
    %plot (abs(new_data));
    %%figure;
    %plot (abs(s_data));
K  = 200;
new_data = fft(s_data);
[~, pos] = sort(new_data(:,i),'descend');
    s_data(pos(1:K),i) = new_data(pos(1:K),i) ;
    figure;
    plot (abs(new_data));
    figure;
    plot (abs(s_data));
actfin_data  = ifft(ifft(new_data));
figure;
plot(abs(actfin_data));
final_data = ifft(ifft(s_data));
figure;
plot(abs(final_data));
 error = norm(final_data - force,2)^2/norm(force,2)^2;
    
    