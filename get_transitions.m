function [EEG] = get_transitions(xpos,ypos,eeg_one,eeg_two,fs)
[searchlog,sinuosity,range,meanspeed] = david_sinuosity3(xpos,ypos,0); %% controltype = 0 for normal, 1 for non-sinuous, 2 for same speed, non-sinuous
[searchlog_control,~,~,~] = david_sinuosity3(xpos,ypos,1);

%% do lfp analysis

%% do coherence during tortuosity/non
%%% tortuosity
sig_tortuous_resamp = resample(double(searchlog),fs,50);
sig_tortuous_resamp_control = resample(double(searchlog_control),fs,50);

for i = 1:size(sig_tortuous_resamp)
    if sig_tortuous_resamp(i,1) == 0
        sig_tortuous_resamp(i,1) = 0;
    else
        sig_tortuous_resamp(i,1) = 1;
    end
end


for i = 1:size(sig_tortuous_resamp_control)
    if sig_tortuous_resamp_control(i,1) == 0
        sig_tortuous_resamp_control(i,1) = 0;
    else
        sig_tortuous_resamp_control(i,1) = 1;
    end
end
%% grab 4s either side centered on when regime changes from non-tort to tort and do
%%% power/coherence
eeg_tortuosity_index = strfind(sig_tortuous_resamp',[0 1]);

for z = 1:size(eeg_tortuosity_index,2)
    if eeg_tortuosity_index(1,z) - (fs*4) < 0 || eeg_tortuosity_index(1,z) + (fs*4) > size(eeg_one,1) %%% check to make sure they're no out of range
        continue
    else
        transition_window_eeg_one(:,z) = eeg_one(eeg_tortuosity_index(1,z)-(fs*4):eeg_tortuosity_index(1,z)+(fs*4),1);
        transition_window_eeg_two(:,z) = eeg_two(eeg_tortuosity_index(1,z)-(fs*4):eeg_tortuosity_index(1,z)+(fs*4),1);
    end
end

EEG.nbchan = 2;
EEG.srate = fs;
EEG.pnts = size(transition_window_eeg_one,1);
EEG.trials = size(transition_window_eeg_one,2);
EEG.data(1,:,:) = transition_window_eeg_one;
EEG.data(2,:,:)= transition_window_eeg_two;


clear transition_window_eeg_one transition_window_eeg_two eeg_tortuosity_index
num_events = length(strfind(sig_tortuous_resamp',[1 0]));

%% CONTROL
% sig_tortuous_resamp_control(sig_tortuous_resamp > 0,:) = 1;

a = 1;
    while a < num_events
    rando = randi([1 size(eeg_one,1)]);
    if sig_tortuous_resamp_control(rando) == 1
        eeg_tortuosity_index(1,a) = rando;
        clear rando
        a = a+1;
    else
    end
    end


for z = 1:size(eeg_tortuosity_index,2)
    if eeg_tortuosity_index(1,z) - (250*4) < 0 | eeg_tortuosity_index(1,z) + (250*4) > size(eeg_one,1)
        continue
    else
        transition_window_eeg_one(:,z) = eeg_one(eeg_tortuosity_index(1,z)-250*4:eeg_tortuosity_index(1,z)+250*4,1);
        transition_window_eeg_two(:,z) = eeg_two(eeg_tortuosity_index(1,z)-250*4:eeg_tortuosity_index(1,z)+250*4,1);
    end
end

EEG.control_data(1,:,:) = transition_window_eeg_one;
EEG.control_data(2,:,:) = transition_window_eeg_two;