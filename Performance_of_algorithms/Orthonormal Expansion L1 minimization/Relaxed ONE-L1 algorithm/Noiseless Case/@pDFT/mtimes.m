function res = mtimes(a,b)

if a.adjoint
    bfull = zeros(a.len,1);
    bfull(a.Omega) = b;
    res = ifft(bfull)*sqrt(a.len);
else
    fb = fft(b)/sqrt(a.len);
    res = fb(a.Omega);
end



    
