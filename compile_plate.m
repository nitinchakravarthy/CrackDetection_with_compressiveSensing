%{
  damagedElementNum=9; %( an array)
  degreeDamage=0.75; %(an array)
%%
clear all;
clc;
method_types = input (['enter \n1: for l1-magic'...
                            '\n2: for Douglas Rachford Method'...
                            '\n3: for SPGL1'...
                            '\n4: for OMP'...
                            '\n5: for CoSaMP'... 
                            '\n6: for gOMP'...
                            '\n7: for Forward Backard Pursuit '...
                            '\n8: for stCSMUSIC '...
                            '\n9: for MUSIC algorithm '...
                            '\n10: for L1 HOMOTOPY'...
                            '\n11: for NNS SP algorithm'...
                            '\n12: for NNS OMP algorithm'...
                            '\n13: for NNS HTP algorithm'...
                            '\n14: for NNS COSAMP algorithm'...
                            '\n15: for all the algorithms\n']);
                            
if method_types == 15
    method_types = 1:14;
end
projectionmatrix_used = input(['In compression use a projection matrix: \n 1 Yes'...
                                                                      '\n 0 NO\n']);
%}
%%
length =8;     %plate length
width =8;       %plate width
seed_length = 1;    %mesh seed length (size of element)
influence_Radius =2;    %depends upon sensor
save influence_Radius.mat influence_Radius
N1=length/seed_length; %along length %aligned with Y-axis
N2=width/seed_length; %along width %aligned with X-axis
[SensList] = CreateSensorList( N2,N1,influence_Radius);

[stiffness, mass, C, freq,~]=generateMindlinPlateData(length,...
                                                        width,seed_length);

[loading] = create_loading_onGrid( freq,2,mass,N1,N2);

A=createFeatureMatrix_plate(SensList,influence_Radius,0);
    damagedElementNums = [ 3, 4, 11, 12];
    degreeDamage = 0.6;
  vDamage=createDamageSignal_plate(damagedElementNums,degreeDamage);
  
  % by the end of the loop each element of the vDamage_group cell should be a (compressed size*dinput_length) size 
  %obs=zeros(size(SensList,2),3);
  for sensor_n = 1: size(SensList,2) 
        
    [ElementListunderInfluence, NodeListunderInfluence]=...
            ElementsUnderInfluenceRadius(SensList(sensor_n),influence_Radius);
        
        %SensLocation_index=find(activeNodeW==SensLocation);
        %Sets of continueous 4 elements 
    DamagedElementSets = Create_DamagedElementSets4(ElementListunderInfluence,...
                                    NodeListunderInfluence);
       z=evaluateOptimumZ(A(:,:,sensor_n),vDamage(:,sensor_n),1);
       figure;plot(1:size(z,1),abs(z));
       [location,severity,l1norm]=identifyDamage_plate(vDamage(:,sensor_n),...
                                          A(:,:,sensor_n),z,DamagedElementSets);
       %{
       for inx = 1:length(location)                               
           location(inx,1)=elementList(location(inx,1));                                                           
           obs(inx,:,sensor_n)=[location(inx,1),severity(inx,1),l1norm(inx,1)]; 
       end
        %}
       obs(sensor_n,:)= [location,severity,l1norm];           
  end
  obs