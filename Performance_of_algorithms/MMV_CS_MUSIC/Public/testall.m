%% Example files to compare several joint sparse recovery algorithm
%% Nov/28/2012 by Jong Chul Ye (jong.ye@kaist.ac.kr), KAIST
%% 

close all

%% Simulation parameter
snr= 30; r_rank = 5;
num_trial = 20;
tau = 1;
n = 128; k=10; 



%% initialization

rand('seed', sum(100*clock));
randn('seed', sum(100*clock));



N_array = [r_rank+10 256]; nN=length(N_array);
m_array = k+1:50; nm = length(m_array);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% perfect ratio

hit_ratio_sa = zeros(nm, nN);
hit_ratio_bs = zeros(nm, nN);
hit_ratio_somp = zeros(nm, nN);
hit_ratio_music = zeros(nm, nN);

% time
time_sa = zeros(nm, nN);
time_bs = zeros(nm, nN);
time_somp = zeros(nm, nN);
time_music = zeros(nm, nN);

%trun = 0;

%% run simulation
for qt = 1 : num_trial      
    tic;
    for N = N_array 
     iN= find(N_array==N);
    for m = m_array   
     im= find(m_array==m);
        if m >= k

           
            % forward
            %[Y,A,X,supp] = forward_fourier(m,n,k,r_rank,snr, N, 0, tau);
  
            [Y,A,X,supp] = forward_gaussian(m,n,k,r_rank,snr, N, 0, tau);
       
            res = Y-A*X; sigma = norm(res(:));
                        
             % CS-MUSIC reconstruction
             t=cputime;
             [actSet, Xh] = CSMUSIC(Y, A, k);  
             time_sa(im,iN) = time_sa(im,iN)+ cputime-t;

             if length(supp)==length(actSet)
                 hit_ratio_sa(im,iN) = hit_ratio_sa(im,iN) +  isempty(setdiff(supp, actSet));
             end

            % SeqCSMUSIC reconstruction 
            t =cputime;
            [actSet, Xh] = SeqCSMUSIC(Y, A, k); 
            time_bs(im,iN) = time_bs(im,iN)+ cputime-t;
            if length(supp)==length(actSet)
                 hit_ratio_bs(im,iN) = hit_ratio_bs(im,iN) + isempty(setdiff(supp, actSet));
            end


            % S-OMP reconstruction 
            t =cputime;
            [actSet, Xh] = SOMP(Y, A, k);
                  time_somp(im,iN) = time_somp(im,iN)+ cputime-t;
            if length(supp)==length(actSet)
                 hit_ratio_somp(im,iN) = hit_ratio_somp(im,iN) + isempty(setdiff(supp, actSet));
            end
             
            
             
             % MUSIC
             t =cputime;
             [actSet, Xh] = MUSIC(Y, A, k);
             time_music(im,iN) = time_music(im,iN)+ cputime-t;
             if length(supp)==length(actSet)
                 hit_ratio_music(im,iN) = hit_ratio_music(im,iN) + isempty(setdiff(supp, actSet));
            end
            
        end             
    end 
    end
    toc;
    
    disp([num2str(qt) ' th simulation']);  
end

hit_ratio_sa = hit_ratio_sa/num_trial;
hit_ratio_bs = hit_ratio_bs/num_trial;
hit_ratio_somp = hit_ratio_somp/num_trial;
hit_ratio_music = hit_ratio_music/num_trial;


time_sa = time_sa/num_trial;
time_bs = time_bs/num_trial;
time_somp = time_somp/num_trial;
time_music = time_music/num_trial;



%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures


plot_chr = {'b-', 'ro:', 'kx-.', 'gs-', 'm-','bx-','r-.', 'ks-'};
method = {'CS-MUSIC', 'SeqCSMUSIC','S-OMP',    'MUSIC'};



for N = N_array  
    iN= find(N_array==N);
    figure;
    hold on
    plot(m_array, squeeze(hit_ratio_sa(:,iN)), plot_chr{1}, 'LineWidth', 1.5);
    plot(m_array, squeeze(hit_ratio_bs(:,iN)), plot_chr{2}, 'LineWidth', 1.5);
    plot(m_array, squeeze(hit_ratio_somp(:,iN)), plot_chr{3}, 'LineWidth', 1.5);
    plot(m_array, squeeze(hit_ratio_music(:,iN)), plot_chr{4}, 'LineWidth', 1.5);
    xlabel('m', 'FontSize', 15);
    ylabel('Success rate', 'FontSize', 15);
    legend(method)
 
end



for N = N_array  
    iN= find(N_array==N);
    figure;
    %  axis([min(m_array), max(m_array), 0 6e-3]);
    hold on
    plot(m_array, squeeze(time_sa(:,iN)), plot_chr{1}, 'LineWidth', 1.5);
    plot(m_array, squeeze(time_bs(:,iN)), plot_chr{2}, 'LineWidth', 1.5);
    plot(m_array, squeeze(time_somp(:,iN)), plot_chr{3}, 'LineWidth', 1.5);
    plot(m_array, squeeze(time_music(:,iN)), plot_chr{4}, 'LineWidth', 1.5);
    xlabel('m', 'FontSize', 15);
    ylabel('Average CPU time (second)', 'FontSize', 15);
    legend(method)
end




%%