function outarray =sin_stretch_event(inarray,lookahead)
% find event start and stretch them out by lookahead samples
% v 1.1 bufix in line 8, == changed to =
outarray=inarray;
% Find the start of the true values
for i= 2:length(inarray)-lookahead
   if and(inarray(i-1)==0, inarray(i)==1)
     outarray(i:i+lookahead)=1;
   end  % if
end % for

