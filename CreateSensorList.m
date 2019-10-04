function [SensList] = CreateSensorList( numElemsX,numElemsY,influence_Radius)
Ycounter=floor(numElemsY/(2*influence_Radius));
Xcounter=floor(numElemsX/(2*influence_Radius));
SensList = zeros(1,(Ycounter)*Xcounter);
SensList_index=1;
for n =0: Ycounter-1 
    for i=1:Xcounter
        SensList(1,SensList_index)= (2*influence_Radius)*n*(numElemsX+1) +...
                                     influence_Radius*(numElemsX+1) + ...
                                     (2*(i-1)*influence_Radius+(influence_Radius+1));
        SensList_index = SensList_index+1;
    end
end

save SensList.mat SensList
end

