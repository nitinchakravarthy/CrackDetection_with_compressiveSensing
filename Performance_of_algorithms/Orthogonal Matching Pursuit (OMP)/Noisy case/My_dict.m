
rand_freq  = zeros(10,10);
for i = 1:10
    x = randperm(100,10);
   rand_freq(i,:) =  x*i;
end

%dict = zeros(1000,100);

for i = 1:10
    x = rand_freq(i,:);
    beam_generator;
   for j = 1:10       
       dict(:,(10*(i-1)+j)) = v_acc(j,:)';
   end       
end