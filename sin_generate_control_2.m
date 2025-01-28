function logout=sin_generate_control_2(login,speedarray,Fs)

% control 2 involves generating new pseudo events placed in less-sinuous
% locations  (allowing for some overlap) and selecting regions with
% similiar speed to actual events (or lower)
% login is logical array with 1 indicating event period
% position input is likely sampled at 50Hz.

printon1=false;               % report internmediate stats for events
printon2=false;               % report output stats for events

overlap=Fs/5;                 % how many samples are allowed to overlap between any login and logout  event
speedbuffer=0.5;              % how many stdev above mean speed are you allowed for control events (lower = more restrictive)
base_event_durn=20;           % duration of base control event in samples
min_duration=10;              % don't let pseudoevent duration go under this many samples
max_duration=Fs*4;            % don't let pseudoevent duration go over this many samples

[numevents,durnevents,btwevents,durnstd,eventstart,eventend]=get_logevent_stats(login);   % get stats on input data

%% get mean speed across all login events

idx=find(login);
speedin=mean(speedarray(idx));              % mean speed during real event
speedin_std=std(speedarray(idx));           % std dev of speed in real event
idx=find(~login);
speedout=mean(speedarray(idx));             % mean speed outside event (just for reference - not used here) 

%% generate pseudo events that are in region not marked as sinuous in login (allowing for overlap)

logout=zeros(length(login),1);                      % initialise logout

durn_this_event=base_event_durn;                    % start with this duration event
i=1;


while i < (length(logout) - max_duration)   % loop through all data
  idx=find(login(i:i+durn_this_event));             % look ahead to see if this region is clear of events
  if length(idx) < overlap                          % if clear of input events (allowing for overlap) 
    rn = normrnd(0,durnstd);                        % generate normal distribution random number
    durn_this_event=round(durn_this_event+rn);      % then add have some random (normal) variation in event length

    if durn_this_event < min_duration; durn_this_event=min_duration;end     % Make sure you don't get too short an event
    if durn_this_event >= max_duration; durn_this_event=max_duration;end     % Make sure you don't get too short an event

    logout(i:i+durn_this_event)=1;                  % set control event here

    i=i+Fs/2;                                       % move ahead to next potential space
  else
    i=i+1;                                          % move ahead 1 space
  end

end
[numevents2,durnevents2,btwevents2,durnstd2,eventstart,eventend]=get_logevent_stats(logout);    % use to describe output and get event start and end

if printon1
  fprintf('#eventsin= %3.0f    #eventsout= %3.0f  \n',numevents,numevents2);
end



%% now remove events where mean speed is higher than speedin + stdev * speedbuffer

for i = 1: length(eventstart)
  if mean(speedarray(eventstart(i):eventend(i))) > (speedin+speedin_std * speedbuffer)           % remove events where speed is too high
    logout(eventstart(i):eventend(i))=0;
  end 
end

idx=find(logout);
controlspeed=mean(speedarray(idx));                 % calculate mean speed during pseudo events 

[numevents3,durnevents3,btwevents3,durnstd3,eventstart,eventend]=get_logevent_stats(logout);   % use to describe output


if printon2
  fprintf('#eventsin= %3.0f    #eventsout= %3.0f speedin= %3.0f  controlspeed= %3.0f \n',numevents,numevents3,speedin,controlspeed);
end

return