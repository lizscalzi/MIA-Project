function logout=sin_generate_control_1(login,Fs)

% control 1 generates new pseudoevents that are the
% average length of the actual events (with variance based on input sd) but are placed in non-sinuous
% locations  (allowing for some overlap)
% login is logical array of actual events with 1 indicating event period
% position input is likely sampled at 50Hz.
% logout is pseudo event array

printon=false;                % report stats for events

overlap=Fs/5;                 % how many samples are allowed to overlap between any login and logout  event
base_event_durn=20;           % duration of base control event in samples
min_duration=10;              % don't let pseudoevent duration go under this many samples
max_duration=Fs*4;            % don't let pseudoevent duration go over this many samples

[numevents,durnevents,btwevents,durnstd,eventstart,eventend]=get_logevent_stats(login);   % get stats on input data


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

if printon
  fprintf('#eventsin= %3.0f    #eventsout= %3.0f  \n',numevents,numevents2);
end
