function vDamage=createMultipleDamageSignal_plate(Damages_multiple )
  % LOADING DATA
    load SensList.mat SensList
    load data.mat data
    load loading.mat loading
    load boundary_conditions.mat boundary_conditions
    load activeDof.mat activeDof
    load set.mat set
    load numactiveDofW.mat numactiveDofW
    t = 0.0:.001:.999;
    numberPoints=1000;
    %load activeNodeW.mat activeNodeW
    %load severityScale.mat severityScale
     numberNodesY   =data(8)+1; %along Y-direction(length)
    numberNodesX   =data(9)+1; %along X-direction(width)
    numberNodes    =numberNodesX*numberNodesY; 
    GDOF          =data(12);
    SensLocation_indices = create_SensLocationindices();
    index_counter=1;
    vDamage=zeros(length(set),length(SensList));
    numberOfDamages = 2;
    DamageQuadrants = [1,4];
    for sensor_n=1:length(SensList)
        SensLocation=SensList(sensor_n);
        SensLocation_index=SensLocation_indices(index_counter,1);   %index in the response
        [stiffness,mass]=generateMultipleDamagedFeature_plate(data,...
                                    boundary_conditions,Damages_multiple);
                                          
            
         [stiffness,mass, damping]=createRayleighClassicalDamping_plate...
                                           (stiffness,mass,activeDof);
        [~, ~, acc_m] = NewmarkMethod(mass,stiffness,damping,loading,t);

        response = acc_m(1:numactiveDofW,:);
        feature_time=response(SensLocation_index,:);
        feature_time=feature_time';
        feature=fft(feature_time,numberPoints);
        feature=awgn(feature,60);
         
        vDamage(:,sensor_n) = feature(set,:);
        index_counter=index_counter+1;
    end
    
end