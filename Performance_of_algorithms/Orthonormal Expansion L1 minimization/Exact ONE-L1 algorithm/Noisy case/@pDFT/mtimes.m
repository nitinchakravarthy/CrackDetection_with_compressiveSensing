function res = mtimes(a,b)
A_dict = a.A_dict;
[N m] = size(A_dict);
Omega  = a.Omega;
if a.adjoint
    
%bfull = zeros(100,1);
 %   bfull(a.Omega) = b;
   
 fb = ifft(A_dict);
    fb = fb(Omega,:);
    res = fb\b;
   
    %res = ifft(bfull)*sqrt(a.len);
    %{
    x =  phi*A_dict;
    x = inv((x')*(x))*(x');
    res = x * b;
    %}    
else
   fb = dftmtx(1000)*A_dict;
   %fb = fb(Omega,:);
 fb = fb * b/sqrt(a.len);
  fb = fb(Omega,:);
    res = fb;
    
    %res = (phi*(a.A_dict)) * b;

end



    
