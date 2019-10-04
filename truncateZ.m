function delZ=truncateZ(z,caseNumber)
  %casesPerElement=10;
  nonzero=caseNumber;
  %nonzero=casesPerElement*(elementNumber-1)+1:1:casesPerElement*(elementNumber);
  zero=setdiff(1:length(z),nonzero);
  delZ=z;
  delZ(zero)=0;
end