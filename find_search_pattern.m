% find_search_pattern
% Load Desiree data as stored in structure built by Rob
% mark regions where animal is moving slower and more tortuous (exploitation)
% Primary output is searchlog (search logic) which is the same length as
% xpos and ypos with exploitation (targeted search) as one, rest as zero.

% Last updated 22/3/24  1.42pm D Bilkey

Fs=100;
lookwindow=1;                       % number of seconds to look ahead
                                    % looking for exploitation =  tortuous, low speed and limited range sections
tortlimit=2;                        % manual setting looking for sections > this
speedlimit=20;                      % manual setting looking for sections < this
rangelimit=50;                      % manual looking for sections < this
autoset=true;                       % ignore above and set thresholds according to threshold
threshold=0.5;                      % stdev from mean for threshold
gapwindow=1;                        % fill gaps between events of less than n seconds
userange=true;                      % use range value to further constrain output

ploton=true;                        % plot paths with events in red
plothits=true;                      % only plot from start of each event. If false plot all blocks
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%

lookahead=lookwindow*Fs;                     % # samples to look ahead for speed
lookahead_low=lookahead/10;                  % # samples to look ahead in resampled location data
Fs_low=Fs/10;                                % Fs for resampled data
pixcm=1;

clear speed
clear tortuosity
clear rangearray
close all

%% load data saved by Rob
load('CTL_Tortuosity_spike_Results.mat');
s=tortuosity_spike_results_ctl;

sz=size(s);
numrecordings=sz(2);
for rnum=1: numrecordings               % recording number  - loop through
    x=s(rnum).xpos;                      % copy x and y data
    y=s(rnum).ypos;
                           
    numsamps=length(x);
    
    %% calculate tortuosity
    
    xds=resample(x,1,10);                               % resample x position down to 1/10 Fs 
    yds=resample(y,1,10);                               % resample y position down to 1/10 Fs 
    
    buffer=zeros(lookahead_low,1);                      % to pad out data at end
    xds=cat(1,xds,buffer);                              % add buffer to end of array
    yds=cat(1,yds,buffer);
    
    for i=1:length(xds)-lookahead_low                           % loop through all position data
      xblock= xds(i:i+lookahead_low);                           % capture block of interest
      yblock= yds(i:i+lookahead_low);
      [sinuousity,range,meanspeed]= getsinuousity2(xblock,yblock,Fs_low);         % get sinuousity measure from this block
      tortuosity(i,1)=sinuousity;                               % returns one sinuousity and one range value for this block
      rangearray(i,1)=range;  
      speed(i,1)=meanspeed;                                     % get mean speed from this block
    end
    
    % resample back to original Fs
    tortuosity=resample(tortuosity,10,1);                       % resample back to same Fs as original arrays
    rangearray=resample(rangearray,10,1);                       % resample back to same Fs as original arrays
    speed=resample(speed,10,1);
    
    
    %% autoset thresholds
    if autoset
      speedlimit=nanmean(speed)-nanstd(speed)*threshold;              % positive threshold decreases speed threshold
      tortlimit=nanmean(tortuosity)+nanstd(tortuosity)*threshold;      % positive threshold increases tortuoisity threshold
      rangelimit=nanmean(rangearray)+nanstd(rangearray)*threshold;    % positive threshold decreases range threshold
    end
    
    %% convert to logical based on thresholds
    speedlog=logical(speed < speedlimit);               % mark slow
    tortlog=logical(tortuosity > tortlimit);            % mark more tortuous
    rangelog=logical(rangearray < rangelimit);          % mark low range
    searchlog=and(tortlog,speedlog);                    % mark exploit
    
    if userange
      searchlog=and(searchlog,rangelog);
    end

    %% stretch every event to cover the lookahead window then fill gaps
    searchlog =sin_stretch_event(searchlog,lookahead);

    searchlog = sin_fill_gaps(searchlog,gapwindow,Fs); % fill gaps less than gapwindow seconds
    
    
    %% graph measures
    f1=figure;
    f1.Position=[50 200 800 600];
    subplot(4,1,1);
    plot(tortuosity);
    yline(tortlimit,'r');
    title('tortuosity');
    subplot(4,1,2);
    plot(speed);
    yline(speedlimit,'r');
    title('speed'); 
    subplot(4,1,3);
    plot(rangearray);
    yline(rangelimit,'r');
    title('range');
    subplot(4,1,4);
    plot(searchlog);
    title('explore');
    
    %% graph paths, marking exploitation in red
    if ploton 
        if plothits
          %% plot all data for hits only
          %hold on
          searchlog_start=sin_find_transitions(searchlog);
          idx=find(searchlog_start);
          numevents=length(idx);
            
          blocksize=1000;                                                   % size of block to plot
          for event=1:numevents           
            f2=figure;
            f2.Position=[800 200 600 600];
            xlim([-50 50]);
            ylim([-50 50]);
            for i = idx(event):idx(event)+blocksize  
                if i < length(searchlog)
                  
                  switch searchlog(i)
                    case 1                                     
                      colorname='r';
                    otherwise
                      colorname='b';
                  end 
                  
                  xvec=x(i:i+1);
                  yvec=y(i:i+1);
                  line(xvec,yvec,'Color',colorname);
                end  % if i
            end
            hold on
            plot(xvec(end),yvec(end),'r*')
            t= ['run = ', num2str(rnum),'   event ', num2str(event), ' of ' num2str(numevents),'     event start = ', num2str(idx(event))];
            title(t);
            pause
            close(f2);
          end  % event
        
        else
            %% plot all data in blocksize sample blocks
            %hold on
            blocksize=1000;
            for block=1:blocksize:numsamps-blocksize                   
              f2=figure;
              f2.Position=[850 200 600 600];
              xlim([-50 50]);
              ylim([-50 50]);
              for i = block:block+blocksize  
                if searchlog(i)==1                                      % or searchlog
                  colorname='r';
                else
                  colorname='b';
                end 
                xvec=x(i:i+1);
                yvec=y(i:i+1);
                line(xvec,yvec,'Color',colorname);
              end
              title('red=exploit');
              pause
              %close all;
            end   % for block
        end   % plothits
    end % ploton
  end % of for loop through recordings


  
