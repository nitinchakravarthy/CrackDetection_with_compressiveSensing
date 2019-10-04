function A=createFeatureMatrix_plate_multipleDamages(SensList,influence_Radius,projectionmatrix_used )

% filename='SensList_plate.xlsx';
%  SensList=xlsread(filename);
  %rewrite the loading statements
  % LOADING DATA
    load data.mat data
    load loading.mat loading
    load boundary_conditions.mat boundary_conditions
    load activeDof.mat activeDof
    load severityScale.mat severityScale
    load numactiveDofW.mat numactiveDofW
    t = 0.0:.001:.999;
    % PARAMETER SETTING
    %numberElements=data(8);
    numberNodesY   =data(8)+1; %along Y-direction(length)
    numberNodesX   =data(9)+1; %along X-direction(width)
    numberNodes    =numberNodesX*numberNodesY; 
    GDOF          =data(12);
    casesPerElement=length(severityScale);
    NumDamagedElementsSets = (2*influence_Radius-1)*(2*influence_Radius-1);
    DamagedElementSets = zeros(NumDamagedElementsSets,4,length(SensList));
    for sensor_n = 1: length(SensList)
    [ElementListunderInfluence, NodeListunderInfluence]=...
            ElementsUnderInfluenceRadius(SensList(sensor_n),influence_Radius);
        
        %SensLocation_index=find(activeNodeW==SensLocation);
        %Sets of continueous 4 elements 
        DamagedElementSets(:,:,sensor_n) = Create_DamagedElementSets4(ElementListunderInfluence,...
                                    NodeListunderInfluence);
    end
    % Number of damages and damage quadrants are predefined(just an assumption)
    numberOfDamages = 2;
    DamageQuadrants = [1,4];
    DamagedElementSets1_full = DamagedElementSets(:,:,DamageQuadrants(1));
    DamagedElementSets1 = DamagedElementSets1_full(1:4,:);
    DamagedElementSets2_full = DamagedElementSets(:,:,DamageQuadrants(2));
    DamagedElementSets2 = DamagedElementSets2_full(1:4,:);
    
    numberTrainingCases=(size(DamagedElementSets1,1)*casesPerElement*...
                            size(DamagedElementSets2,1)*casesPerElement); 
    numberPoints=1000; % parameters regarding fourier transform
    A=zeros(numberPoints,numberTrainingCases,length(SensList));
    
   h = waitbar(0,'Please wait while creating Dictionary (multiple Damages)...');
    
    SensLocation_indices = create_SensLocationindices();
    index_counter=1;
    for layer=1:length(SensList)
        SensLocation=SensList(layer);  
        SensLocation_index=SensLocation_indices(index_counter,1);   %index in the response
        
        % EVALUATE FEATURE MATRIX
    
        for numElement1=1:size(DamagedElementSets1,1)
            for severityCount1=1:casesPerElement
               for numElement2 = 1: size(DamagedElementSets2,1)
                   for severityCount2 = 1:casesPerElement
                        %Waitbar
                       waitbar_length = length(SensList)...
                                        *size(DamagedElementSets1,1)...
                                        *casesPerElement...
                                        *size(DamagedElementSets2,1)...
                                        *casesPerElement;
                        waitbar_filled = casesPerElement*size(DamagedElementSets2,1)...
                                        *casesPerElement*(numElement1-1)...
                                        +casesPerElement*size(DamagedElementSets2,1)...
                                        *severityCount1...
                                        +casesPerElement*(numElement2-1)...
                                        +severityCount2 ;
                        waitbar(waitbar_filled/waitbar_length,h);
                
                
                        % Evaluate Features
                        
                % For a plate single element damage does not make a lot of 
                % difference physically. So a set of consequitve four 
                % elements are taken as damaged elements.
                        
                        Damages_multiple = {DamagedElementSets1(numElement1,:),...
                                            severityScale(severityCount1);...
                                            DamagedElementSets2(numElement2,:),...
                                            severityScale(severityCount2)};
                        [stiffness, mass]=generateMultipleDamagedFeature_plate(...
                                            data,boundary_conditions,Damages_multiple);
                
                        [stiffness,mass, damping]=createRayleighClassicalDamping_plate...
                                                    (stiffness,mass,activeDof);
                
                        [~, ~, acc_m] = NewmarkMethod(mass,stiffness,damping,loading,t);
                        response = acc_m(1:numactiveDofW,:);
                
                        feature_time=response(SensLocation_index,:);
                        feature_time=feature_time';
                        feature=fft(feature_time,numberPoints);
                        feature=awgn(feature,60);
                        Aindex=((casesPerElement*(numElement1-1)+severityCount1)-1)...
                                *(size(DamagedElementSets2,1)*casesPerElement)+...
                                (casesPerElement*(numElement2-1)+severityCount2);
                        A(:,Aindex,layer)=feature;
                   end    
               end
            end
        end
        index_counter=index_counter+1;
    end
    if ~projectionmatrix_used 
     %   use [A,set]=reduceSignal(A,M); where M is the size of compression
     %Default size =10
    [A,set]=reduceSignal(A,10);
    save set.mat set
    else
     %use [A,set]=reduceSignal(A,M); where M is the size of compression
     %Default size =10
    [A,set]=reduceSignal_projmatrixused(A);
    save set.mat set;
    end
    close(h); 
end

