function A=createFeatureMatrix_plate(SensList,influence_Radius,projectionmatrix_used )

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
    [ElementListunderInfluence, NodeListunderInfluence]=...
            ElementsUnderInfluenceRadius(SensList(1),influence_Radius);
        
        %SensLocation_index=find(activeNodeW==SensLocation);
        %Sets of continueous 4 elements 
        DamagedElementSets = Create_DamagedElementSets4(ElementListunderInfluence,...
                                    NodeListunderInfluence);
    numberTrainingCases=(size(DamagedElementSets ,1)*casesPerElement); 
    numberPoints=1000; % parameters regarding fourier transform
    A=zeros(numberPoints,numberTrainingCases,length(SensList));
    
   h = waitbar(0,'Please wait while creating training data...');
    
    SensLocation_indices = create_SensLocationindices();
    index_counter=1;
    for layer=1:length(SensList)
        SensLocation=SensList(layer);  
        SensLocation_index=SensLocation_indices(index_counter,1);   %index in the response
        % Elements and nodes in the influence in ascending order
        [ElementListunderInfluence, NodeListunderInfluence]=...
            ElementsUnderInfluenceRadius(SensLocation,influence_Radius);
        
        %SensLocation_index=find(activeNodeW==SensLocation);
        %Sets of continueous 4 elements 
        DamagedElementSets = Create_DamagedElementSets4(ElementListunderInfluence,...
                                    NodeListunderInfluence);
        % EVALUATE FEATURE MATRIX
    
        for numElement=1:length(DamagedElementSets)
            for severityCount=1:casesPerElement
               
                %Waitbar
                       
                        waitbar_length = length(SensList)...
                                        *length(DamagedElementSets)...
                                        *casesPerElement;
                                        
                        waitbar_filled = casesPerElement*length(DamagedElementSets)...
                                        *(layer-1)...
                                        +casesPerElement*(numElement-1)...
                                        +severityCount ;
                        waitbar(waitbar_filled/waitbar_length,h);
                
                
                        % Evaluate Features
                        
                % For a plate single element damage does not make a lot of 
                % difference physically. So a set of consequitve four 
                % elements are taken as damaged elements.
                
                [stiffness, mass]=generateDamagedFeature_plate(data,boundary_conditions,...
                                           DamagedElementSets(numElement,:),...
                                           severityScale(severityCount));
                
               [stiffness,mass, damping]=createRayleighClassicalDamping_plate...
                                           (stiffness,mass,activeDof);
                
                [~, ~, acc_m] = NewmarkMethod(mass,stiffness,damping,loading,t);
                response = acc_m(1:numactiveDofW,:);
                
                feature_time=response(SensLocation_index,:);
                feature_time=feature_time';
                feature=fft(feature_time,numberPoints);
               % feature=awgn(feature,60);
                Aindex=casesPerElement*(numElement-1)+severityCount;
                A(:,Aindex,layer)=feature;
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

