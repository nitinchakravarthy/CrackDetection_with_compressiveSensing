     '
     'Ax - b is small enough, given atol, btol                  '
     'The least-squares solution is good enough, given atol     '
     'The estimate of cond(Abar) has exceeded conlim            '
     'Ax - b is small enough for this machine                   '
     'The least-squares solution is good enough for this machine'
     'Cond(Abar) seems to be too large for this machine         '
     'The iteration limit has been reached                      '];

wantvar= nargout >= 6;
if wantvar, var = zeros(n,1); end

if show
   disp(' ')
   disp('LSQR            Least-squares solution of  Ax = b')
   str1 = sprintf('The matrix A has %8g rows  and %8g cols', m, n);
   str2 = sprintf('damp = %20.14e    wantvar = %8g', damp,wantvar);
   str3 = sprintf('atol = %8.2e                 conlim = %8.2e', atol, conlim);
   str4 = sprintf('btol = %8.2e                 itnlim = %8g'  , btol, itnlim);
   disp(str1);   disp(str2);   disp(str3);   disp(str4);
end

itn    = 0;		 istop  = 0;        nstop  = 0;
ctol   = 0;		 if conlim > 0, ctol = 1/conlim; end;
anorm  = 0;		 acond  = 0;
dampsq = damp^2; ddnorm = 0;        res2   = 0;
xnorm  = 0;	     xxnorm = 0;        z      = 0;
cs2    = -1;     sn2    = 0;

% Set up the first vectors u and v for the bidiagonalization.

% These satisfy  beta*u = b,  alfa*v = A'u.

u      = b(1:m);	x    = zeros(n,1);
alfa   = 0;		beta = norm( u );
if beta > 0
   u = (1/b