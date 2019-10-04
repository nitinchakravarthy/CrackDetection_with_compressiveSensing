function [locations1,severity1,locations2,severity2,reliability] = ...
         identifyMultipleDamages_plate(vDamage,A,z,DamagedElementSets1,...
                                        DamagedElementSets2)

% A              =  corresponding Dictionary
% vDamage        =  Signal responses of a particular sensor for different 
%                   damage locations and magnitudes 
% z              =  Optimum Z for all the vectors in vDamage_group

numDamages =2;

load severityScale.mat severityScale
casesPerDamage=size(severityScale,2);
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
figure;
plot(safetyCoefficient);
%l1norm=l1norm/max(safetyCoefficient);

%% EXTRACTING DAMAGE PROPERTIES

numberCasesPerEach1stDamage = casesPerDamage*size(DamagedElementSets2,1);
numberCasesPerEach1stDamageLoc = casesPerDamage*numberCasesPerEach1stDamage;
 
FirstLocation       = ceil(damageCase/numberCasesPerEach1stDamageLoc); 
locations1          = DamagedElementSets1(FirstLocation,:);
 
CasesAfterLoc1      = damageCase-(FirstLocation-1)*numberCasesPerEach1stDamageLoc;
FirstSeverity       = ceil(CasesAfterLoc1/numberCasesPerEach1stDamage); 
severity1           = severityScale(FirstSeverity);
 
CasesAfterLoc1Sev1  = damageCase-(FirstLocation-1)*numberCasesPerEach1stDamageLoc...
                        -(FirstSeverity-1)*numberCasesPerEach1stDamage;
SecondLocation      = ceil(CasesAfterLoc1Sev1/casesPerDamage);
locations2          = DamagedElementSets2(SecondLocation,:);
 
SecondSeverity      = damageCase-(FirstLocation-1)*numberCasesPerEach1stDamageLoc...
                        -(FirstSeverity-1)*numberCasesPerEach1stDamage...
                        -(SecondLocation-1)*casesPerDamage;
severity2           = severityScale(SecondSeverity);