function [redSignal,set]=reduceSignal(orSignal,numPoints)
  if nargin<2
     numPoints=10;
  end   
  nCases=size(orSignal,2);
  mag=abs(orSignal);
  set=[];
  
  for layer=1:size(orSignal,3)
      for count=1:nCases
          [~,p_indx]=sort(mag(:,count,layer),'descend');
          points=p_indx(1:numPoints);
          set=union(set,points);
      end
  end
  redSignal=orSignal(set,:,:);
end