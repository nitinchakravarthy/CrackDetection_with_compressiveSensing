function [location,loadingNumber,unloadingNumber,loadMag,reliability]=identifyImpact_plate(vDamage,A,z,l_num)

% A              =  corresponding Dictionary
% vDamage        =  Signal responses of a particular sensor for different 
%                   damage locations and magnitudes 
% z              =  Optimum Z for all the vectors in vDamage_group
    load L_w.mat L_w
    load UNL_w.mat UNL_w
    load activeDof.mat activeDof

casesPerElement=size(L_w,2)*size(UNL_w,2);
numberTrainingCases=size(A,2);
loadingNum  = length(L_w);
unloadingNum = length(UNL_w);                                

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
[sorted_z,pos] = sort(z,'descend');
damageCase = zeros(1,l_num);
loadMag =  zeros(1,l_num);
frequency_scale = zeros(1,l_num);
location = zeros(1,l_num);
for i = 1:1:l_num
damageCase(i) = pos(i);
end

for i = 1:1:l_num
location(i)=ceil(damageCase(i)/casesPerElement);

casesAfterloc = damageCase(i)-casesPerElement*(location(i)-1);
loadingNumber=ceil(casesAfterloc/unloadingNum);

unloadingNumber = damageCase(i)-casesPerElement*(location(i)-1)-unloadingNum*(loadingNumber-1);
loadMag(i) = abs(max((z(damageCase(i)))));
location(i)  = activeDof(location(i));
end

%{
%% EXTRACTING DAMAGE PROPERTIES

numberCasesPerEach1stDamage = casesPerDamage*length(elementList_Damage2);
numberCasesPerEach1stDamageLoc = casesPerDamage*numberCasesPerEach1stDamage;
 
location1           = ceil(damageCase/numberCasesPerEach1stDamageLoc); 
FirstLocation       = elementList_Damage1(location1);
 
CasesAfterLoc1      = damageCase-(location1-1)*numberCasesPerEach1stDamageLoc;
severity1           = ceil(CasesAfterLoc1/numberCasesPerEach1stDamage); 
FirstSeverity       = severityScale(severity1);
 
CasesAfterLoc1Sev1  = damageCase-(location1-1)*numberCasesPerEach1stDamageLoc...
                        -(severity1-1)*numberCasesPerEach1stDamage;
location2           = ceil(CasesAfterLoc1Sev1/casesPerDamage);
SecondLocation      = elementList_Damage2(location2);
 
severity2           = damageCase-(location1-1)*numberCasesPerEach1stDamageLoc...
                        -(severity1-1)*numberCasesPerEach1stDamage...
                        -(location2-1)*casesPerDamage;
SecondSeverity      = severityScale(severity2);
%}
end