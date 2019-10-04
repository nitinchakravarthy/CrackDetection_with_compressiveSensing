function [D,VV,VVV]=...
          calcEigenvectors...
         (mass,stiffness,activeDof,activeNodeW,numberOfModes,data)
  [V,D] = eig(stiffness(activeDof,activeDof),...
                   mass(activeDof,activeDof)); 
  D = diag(sqrt(D))*sqrt(data(1)*data(2))*sqrt(data(7)/(data(5)/2.6)); %for Mindlin
  [D,ii] = sort(D); ii = ii(1:numberOfModes); 
  VV = V(:,ii);
  numberNodes=(data(1)/+1)*(data(2)+1);
  VVV(1:numberNodes,1:numberOfModes)=0;
    for i=1:numberOfModes
        VVV(activeNodeW,i)=VV(1:size(activeNodeW),i);
    end
end