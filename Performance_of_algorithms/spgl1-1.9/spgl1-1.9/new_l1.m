
    % -----------------------------------------------------------
    % Solve the basis pursuit (BP) problem in COMPLEX variables:
    %
    %    minimize ||z||_1 subject to Az = b
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the basis pursuit (BP) problem in COMPLEX variables:\n');
    fprintf('%%                                                    \n');
    fprintf('%%   minimize ||z||_1 subject to Az = b               \n');
    fprintf('%%                                                    \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

  

    % Create partial Fourier operator with rows idx
    idx = randperm(n); idx = idx(1:m);
    opA = @(x,mode) partialFourier(idx,n,x,mode);

    % Create sparse coefficients and b = 'A' * z0;
    z0 = zeros(n,1);
    z0(p) = randn(k,1) + sqrt(-1) * randn(k,1);
    b = opA(z0,1);
    
    opts = spgSetParms('verbosity',1);
    z = spg_bp(opA,b,opts);
    
    figure(1); subplot(2,4,3);
    plot(1:n,real(z),'b+',1:n,real(z0),'bo', ...
         1:n,imag(z),'r+',1:n,imag(z0),'ro');
    legend('Recovered (real)', 'Original (real)', ...
           'Recovered (imag)', 'Original (imag)');
    title('(c) Complex Basis Pursuit');
    
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(c).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');
    