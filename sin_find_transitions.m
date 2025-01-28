function outarray=sin_find_transitions(inarray)
% mark start regions of logical on
arraylen=length(inarray);
outarray=zeros(arraylen,1);
for i=1:arraylen-1
  if and(inarray(i)==0,inarray(i+1)==1)
     outarray(i+1)=1;
  end
end
