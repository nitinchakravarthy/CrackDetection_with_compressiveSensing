 A = zeros(100);
 for i =1:100 
     sum = zeros(1,100);
     for j = 1:10 
         
         sum = sum + final_acc((i-1)*10+j,:);
     end
     A(i,:) = sum; 
 end