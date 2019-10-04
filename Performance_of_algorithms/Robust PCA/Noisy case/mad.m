function [MAD]= mad(Z,dim)

%Median Absolute Deviation
%if the dimension is not specified or dim =1 MAD is a row vector 
%containing the MAD of each column
%dim = 2 gives a colum vector containing the MAD of each row
if nargin<2 
   dim = 1;
end
if dim==1
   MAD=median(abs(Z - repmat(median(Z,dim),size(Z,dim),1)),dim)*1.4826; 
else
   MAD=median(abs(Z - repmat(median(Z,dim),1,size(Z,dim))),dim)*1.4826; 
end












