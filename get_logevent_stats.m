function [numevents,durnevents,btwevents,durnstd,eventstart,eventend]=get_logevent_stats(inlog)
% returns number of events, duratino of events (samples), samples between
% events, stdev of duration (samples) event start and event end indices
len=length(inlog);
inlog(1)=0;
inlog(end)=0;

logout(1)=0;
logout(end)=0;
transition=diff(inlog);                     % find event starts and ends
eventstart=find(transition >0);             % idx to start of events
eventend= find(transition < 0);             % idx to end of events

if and(length(eventstart)>0,length(eventend)>0)
    numevents=length(eventstart);
    idx=find(inlog);
    num_on=length(idx);
    durnevents=round(num_on/numevents);    % average duration of events in samples
    num_off=len-num_on;
    btwevents=round(num_off/numevents);
    
    for i=1:numevents
      durn(i)=eventend(i)-eventstart(i);
    end
    durnstd=std(durn);
else
    numevents=0;
    durnevents=0;
    btwevents=0;
    durnstd=0;
end

return