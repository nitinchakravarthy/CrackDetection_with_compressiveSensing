function  res = pDFT(len, Omega)

if nargin < 2
	Omega = 1:len;
end

res.adjoint = 0;
res.Omega = Omega;
res.len = len;

res = class(res,'pDFT');

end
