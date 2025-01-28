
function [searchlog,sinuousity,range,meanspeed] = david_sinuosity3(xpos,ypos,controltype)
%%% Takes position data and returns:
%%% searchlog: binary "search" and "exploit" sections
%%% sinuosity: (a VTI-like measure)
%%% range: The space covered in a two second window episode
%%% meanspeed: the mean running speed
%%% Version 2 is modified to  calculate range based on eccentricity of path
%%% Version 3 adds control options by detecting various other sections of
%%% searchlog is based on range and speed.

%%% assumes Fs of input is 50Hz
%%% D Bilkey 22/5/24
%%% Modified to take controltype as an argument rather than a hard-code, R
%%% Munn 10/6/24

% controltype=0;                      % 0 = run normally  (control is off)
                                    % 1 = detects non-sinuous regions
                                    % 2 = detects events with similar speed in non
                                    %     sinuous region
                               

possample=50;                       % position sampling rate
lookwindow=2;                       % number of seconds to look ahead
resamplefactor=5;                   % downsample data by n to get 10 samples per second
                                    % looking for exploitation =  tortuous, low speed and limited range sections
tortlimit=2;                        % manual setting looking for sections > this
speedlimit=20;                      % manual setting looking for sections < this
rangelimit=25;                      % manual looking for sections < this
autoset=false;                      % ignore above and set thresholds according to threshold
threshold=0.5;                      % stdev from mean for threshold
gapwindow=1;                        % fill gaps between events of less than n seconds
userange=true;                      % use range value to further constrain output

ploton=true;                        % plot paths with events in red
plothits=true;                      % only plot from start of each event. If false plot all blocks
Fs = possample;              
%%%%%%%%%%%%%%%%%%%%%%%%%%%

lookahead=lookwindow*possample;                                 % # samples to look ahead for speed
lookahead_low=lookahead/resamplefactor;                         % # samples to look ahead in resampled location data
Fs_low=possample/resamplefactor;                                % Fs for resampled data
pixcm=1;

clear speed
clear tortuosity
clear rangearray
clear speedarray
clear tortuosityarray
close all


x = xpos;
y = ypos;
%% calculate tortuosity, range, speed

xds=resample(x,1,resamplefactor);                               % resample x position down to 1/10 Fs 
yds=resample(y,1,resamplefactor);                               % resample y position down to 1/10 Fs 

buffer=zeros(lookahead_low,1);                                  % to pad out data at end
xds=cat(1,xds,buffer);                                          % add buffer to end of array
yds=cat(1,yds,buffer);

for i=1:length(xds)-lookahead_low                               % loop through all position data
  xblock= xds(i:i+lookahead_low);                               % capture block of interest
  yblock= yds(i:i+lookahead_low);
  [sinuousity,range,meanspeed]= getsinuousity6(xblock,yblock,Fs_low);         % get sinuousity etc. measure from this block
  tortuosityarray(i,1)=sinuousity;                              % returns one sinuousity and one range value for this block
  rangearray(i,1)=range;  
  speedarray(i,1)=meanspeed;                                    % get mean speed from this block
end  % for

% resample back to original Fs
tortuosityarray=resample(tortuosityarray,resamplefactor,1);                       % resample back to same Fs as original arrays
rangearray=resample(rangearray,resamplefactor,1);                       % resample back to same Fs as original arrays
speedarray=resample(speedarray,resamplefactor,1);


%% autoset thresholds
if autoset
  speedlimit=nanmean(speedarray)-nanstd(speedarray)*threshold;                  % positive threshold decreases speed threshold
  tortlimit=nanmean(tortuosityarray)+nanstd(tortuosityarray)*threshold;         % positive threshold increases tortuoisity threshold
  rangelimit=nanmean(rangearray)+nanstd(rangearray)*threshold;                  % positive threshold decreases range threshold
end  % autoset

%% convert to logical based on thresholds
speedlog=logical(speedarray < speedlimit);               % mark slow
tortlog=logical(tortuosityarray > tortlimit);            % mark more tortuous
rangelog=logical(rangearray < rangelimit);               % mark low range

%% stretch every event to cover the lookahead window then fill gaps
speedlog =sin_stretch_event(speedlog,lookahead);
rangelog =sin_stretch_event(rangelog,lookahead);
tortlog =sin_stretch_event(tortlog,lookahead);

searchlog=and(rangelog,speedlog);                    % mark exploit using range and speed

%% return mean speed range and sinuosity values
sinuosity=nanmean(tortuosityarray);
range=nanmean(rangearray);
meanspeed=nanmean(speedarray);


%searchlog =sin_stretch_event(searchlog,lookahead);

searchlog = sin_fill_gaps(searchlog,gapwindow,Fs); % fill gaps less than gapwindow seconds


%% if controltype > 0 then replace searchlog with a new one that marks control events

switch controltype
  case 0
    % do nothing
  case 1
    searchlog=sin_generate_control_1(searchlog,Fs);   % move events to non-sinuous regions
  case 2
    searchlog=sin_generate_control_2(searchlog,speedarray,Fs); % move events to non-sinuous regions with simiar speed
  case 3

 end
   


end


