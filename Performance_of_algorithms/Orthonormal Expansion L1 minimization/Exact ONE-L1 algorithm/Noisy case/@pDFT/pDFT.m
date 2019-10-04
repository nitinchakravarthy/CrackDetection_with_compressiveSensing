 function  res = pDFT(len, Omega,phi,A_dict)

if nargin < 2
	Omega = 1:len;
end

res.adjoint = 0;
res.Omega = Omega;
res.len = len;
res.phi = phi;
res.A_dict = A_dict;

res = class(res,'pDFT');

end
