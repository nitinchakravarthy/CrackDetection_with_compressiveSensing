% Guyan's  Static Reduction
% to get the Mass and Stiffness matrices of the structural model
function [Mr, Kr, Cr, Pr] = guyan(M, K, C, P,pts,x,y)
% Re-arranging the stiffness and mass matrices into four sections
% the mm (master-master), ss (slave-slave), sm (slave-master) and ms
% (master-slave)
if length(y) >1
    for t=1:pts
        for i=1:length(x)
            Kmm(i,:)=K(x(i),x);
            Mmm(i,:)=M(x(i),x);
            Cmm(i,:)=C(x(i),x);
            Pmm(i,t)=P(x(i),t);
        end
        for i=1:length(y)
            Kss(i,:)=K(y(i),y);
            Mss(i,:)=M(y(i),y);
            Css(i,:)=C(y(i),y);
            Pss(i,t)=P(y(i),t);
        end
        for i=1:length(y)
            Ksm(i,:)=K(y(i),x);
            Msm(i,:)=M(y(i),x);
            Csm(i,:)=C(y(i),x);
        end
        
        Kms = Ksm';Mms = Msm';Cms = Csm';
        % transfer matrix that will reduce the matrices
        T = [eye(length(x)) ; -inv(Kss) * Ksm];
        % Placing the matric sections in the appropriate spots.
        Ma = [Mmm Mms ; Msm Mss];
        Ka = [Kmm Kms ; Ksm Kss];
        Ca = [Cmm Cms ; Csm Css];
        Pa = [Pmm(:,t) ; Pss(:,t)];
        % Calculating the reduced mass and stiffness matrices.
        Mr = T' * Ma * T;
        Kr = T' * Ka * T;
        Cr = T' * Ca * T;
        Pr(:,t) = T' * Pa ;
        
    end
else Mr =M; Kr =K; Cr=C; Pr=P;
end