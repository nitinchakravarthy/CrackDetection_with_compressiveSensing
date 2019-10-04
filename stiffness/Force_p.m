function F = Force_p(loc,force_data,p)

GDof = force_data.d1;
numberElements=force_data.d2;
elementNodes=force_data.d3;
numberNodes=force_data.d4;
nodeCoordinates=force_data.d5;
activeDof=force_data.d6;
% load
% p = -1*sin(1.3*t);
P=zeros(numberElements,1);
% P(23)=1; P(28)=1; P(73)=1;P(78)=1;
P(loc)=ones(length(loc),1);
P=P*p;

for i=1:length(p)
Force(:,i)=...
formForceVectorMindlinQ4(GDof,numberElements,...
elementNodes,numberNodes,nodeCoordinates,P(:,i));
F(:,i)=Force(activeDof,i);
end
