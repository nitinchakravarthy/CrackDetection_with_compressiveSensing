function [location,severity,reliability]=identifyDamage(vDamage,A,z)

% A              =  corresponding Dictionary
% vDamage        =  Signal responses of a particular sensor for different 
%                   damage locations and magnitudes 
% z              =  Optimum Z for all the vectors in vDamage_group

load severityScale.mat severityScale
casesPerElement=size(severityScale,2);
numberTrainingCases=size(A,2);

%% EVALUATING FEATURE WITH LOWEST SAFTEY COEFFICIENT
safetyCoefficient=zeros(numberTrainingCases,1);
for caseNumber=1:numberTrainingCases
    %elementNumber=ceil(caseNumber/casesPerElement);
    rMean=norm(vDamage-A*truncateZ(z,caseNumber));
    %rSum=0;
    %for count=1:q
        %phi=randn(size(A,1));
        %Amodified=phi*A;
        %Amodified=A;
        %vModified=phi*vDamage;
        %vModified=vDamage;
        %rSum=rSum+norm(vModified-Amodified*truncateZ(z,caseNumber));
    %end
    %rMean=rSum/q;
    safetyCoefficient(caseNumber)=rMean;
end
[l1norm,damageCase]=min(safetyCoefficient);



others=setdiff(1:length(safetyCoefficient),damageCase);
%reliability=(mean(safetyCoefficient)-l1norm)/mean(safetyCoefficient);
%reliability=std(safetyCoefficient(others))/std(safetyCoefficient);
reliability=mean(safetyCoefficient(others))/std(safetyCoefficient(others));
%figure;
%plot(safetyCoefficient);
%l1norm=l1norm/max(safetyCoefficient);

%% EXTRACTING DAMAGE PROPERTIES

location        =ceil(damageCase/casesPerElement);
severityNumber  =damageCase-casesPerElement*(location-1);
severity        =severityScale(severityNumber);

end