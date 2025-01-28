% test_david_sinuosity2
% tests david_sinuosity2 as constructed by Rob and monified to v2 by David
% Load Desiree data as stored in structure built by Rob
% mark regions where animal is moving slower and more tortuous (exploitation)
% Primary output is searchlog (search logic) which is the same length as
% xpos and ypos with exploitation (targeted search) as one, rest as zero.
% looking for exploitation =  tortuous, low speed and limited range sections

% Last updated 21/4/24  D Bilkey


ploton=false;                            % plot paths with events in red
plothits=false;                         % only plot from start of each event. If false plot all blocks
Fs=50;                                  % sampling rate for position data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%% load data saved by Rob
load('CTL_Tortuosity_spike_Results.mat');
s=tortuosity_spike_results_ctl;

sz=size(s);
numrecordings=sz(2);
for rnum=1: numrecordings               % recording number  - loop through
    xpos=s(rnum).xpos;                      % copy x and y data
    ypos=s(rnum).ypos;
               
    numsamps=length(xpos);

   [searchlog,sinuousity,range,meanspeed] = david_sinuosity2(xpos,ypos);
    
    
   
    
    %% graph paths, marking exploitation in red
    if ploton 
        if plothits
          %% plot all data for hits only
          %hold on
          searchlog_start=sin_find_transitions(searchlog);
          idx=find(searchlog_start);
          numevents=length(idx);
            
          blocksize=500;                                                   % size of block to plot
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
                  
                  xvec=xpos(i:i+1);
                  yvec=ypos(i:i+1);
                  line(xvec,yvec,'Color',colorname);
                end  % if i
            end
            hold on
            plot(xvec(end),yvec(end),'r*')
            t= ['run = ', num2str(rnum),'   event ', num2str(event), ' of ' num2str(numevents),'     event start = ', num2str(idx(event))];
            title(t);
            %pause
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
                xvec=xpos(i:i+1);
                yvec=ypos(i:i+1);
                line(xvec,yvec,'Color',colorname);
                hold on
                if mod(i,Fs)==0
                  plot(xpos(i),ypos(i),'k*');
                end
              end
              hold on
              plot(xvec(end),yvec(end),'r*')
             
              t= ['red=exploit  red asterisk = end     black = 1 second interval     start block= ', num2str(block)];
              title(t);
              %pause
              keyin = input('press x to exit return to continue','s');
              if keyin=='x'
                error('terminated');
              end
              close (f2);
            end   % for block

        end   % plothits
    end % ploton
  end % of for loop through recordings


