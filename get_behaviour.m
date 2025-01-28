function [results] = get_behaviour(xpos,ypos,eeg_one,eeg_two,fs,id)
[searchlog,sinuosity,range,meanspeed] = david_sinuosity3(xpos,ypos,0); %% controltype = 0 for normal, 1 for non-sinuous, 2 for same speed, non-sinuous
[searchlog_control,~,~,meanspeed_control] = david_sinuosity3(xpos,ypos,1);

%% LOW TO HIGH
fs = 50;
%%% grab 4s either side centered on when regime changes from non-tort to tort and do
%%% power/coherence
pos_index_low_to_high = strfind(searchlog',[0 1]);

for z = 1:size(pos_index_low_to_high,2)
    if pos_index_low_to_high(1,z) - (50*4) <= 0 || pos_index_low_to_high(1,z) + (50*4) > size(xpos,1) %%% check to make sure they're no out of range
        continue
    else
        posx_window = xpos(pos_index_low_to_high(1,z)-(fs*4):pos_index_low_to_high(1,z)+(fs*4),1);
        posy_window = ypos(pos_index_low_to_high(1,z)-(fs*4):pos_index_low_to_high(1,z)+(fs*4),1);
        post_window = linspace(-4,4,size(posx_window,1))';
        speed_low_to_high_transition_pre(z,:) = speed2D(posx_window(1:size(posx_window,1)/2,1),posy_window(1:size(posx_window,1)/2,1),post_window(1:size(posx_window,1)/2,1));
        speed_low_to_high_transition_post(z,:) = speed2D(posx_window(size(posx_window,1)/2:end,1),posy_window(size(posx_window,1)/2:end,1),post_window(size(posx_window,1)/2:end,1));
    end
end

for z = 1:size(speed_low_to_high_transition_pre,1)
    if sum(speed_low_to_high_transition_pre(z,:)) == 0
        preremove(z) = false;
    else
        preremove(z) = true;
    end
end

for z = 1:size(speed_low_to_high_transition_post,1)
    if sum(speed_low_to_high_transition_post(z,:)) == 0
        postremove(z) = false;
    else
        postremove(z) = true;
    end
end
results.animal = id;
results.speed_low_to_high_transition_post = speed_low_to_high_transition_post(postremove,:);
results.speed_low_to_high_transition_pre= speed_low_to_high_transition_pre(preremove,:);

results.low_high_pretransition_mean = mean(results.speed_low_to_high_transition_pre,2);
results.low_high_postransition_mean = mean(results.speed_low_to_high_transition_post,2);

clear posx_window posy_window post_window preremove postremove
%%% CONTROL
pos_index_low_to_high_control = strfind(searchlog_control',[0 1]);


for z = 1:size(pos_index_low_to_high_control,2)
    try
    if pos_index_low_to_high_control(1,z) - (50*4) <= 0 || pos_index_low_to_high_control(1,z) + (50*4) > size(xpos,1) %%% check to make sure they're no out of range
        continue
    else
        posx_window = xpos(pos_index_low_to_high_control(1,z)-(fs*4):pos_index_low_to_high_control(1,z)+(fs*4),1);
        posy_window = ypos(pos_index_low_to_high_control(1,z)-(fs*4):pos_index_low_to_high_control(1,z)+(fs*4),1);
        post_window = linspace(-4,4,size(posx_window,1))';
        speed_low_to_high_transition_pre_control(z,:) = speed2D(posx_window(1:size(posx_window,1)/2,1),posy_window(1:size(posx_window,1)/2,1),post_window(1:size(posx_window,1)/2,1));
        speed_low_to_high_transition_post_control(z,:) = speed2D(posx_window(size(posx_window,1)/2:end,1),posy_window(size(posx_window,1)/2:end,1),post_window(size(posx_window,1)/2:end,1));
    end

catch
    keyboard
end
end

for z = 1:size(speed_low_to_high_transition_pre_control,1)
    if sum(speed_low_to_high_transition_pre_control(z,:)) == 0
        preremove(z) = false;
    else
        preremove(z) = true;
    end
end

for z = 1:size(speed_low_to_high_transition_post_control,1)
    if sum(speed_low_to_high_transition_post_control(z,:)) == 0
        postremove(z) = false;
    else
        postremove(z) = true;
    end
end

results.speed_low_to_high_transition_post_control = speed_low_to_high_transition_post_control(postremove,:);
results.speed_low_to_high_transition_pre_control = speed_low_to_high_transition_pre_control(preremove,:);

results.low_high_pretransition_mean_control = mean(results.speed_low_to_high_transition_pre_control,2);
results.low_high_postransition_mean_control = mean(results.speed_low_to_high_transition_post_control,2);
clear posx_window posy_window post_window preremove postremove
%% HIGH TO LOW

pos_index_high_to_low = strfind(searchlog',[1 0]);

for z = 1:size(pos_index_high_to_low,2)
    if pos_index_high_to_low(1,z) - (50*4) <= 0 || pos_index_high_to_low(1,z) + (50*4) > size(xpos,1) %%% check to make sure they're no out of range
        continue
    else
        posx_window = xpos(pos_index_high_to_low(1,z)-(fs*4):pos_index_high_to_low(1,z)+(fs*4),1);
        posy_window = ypos(pos_index_high_to_low(1,z)-(fs*4):pos_index_high_to_low(1,z)+(fs*4),1);
        post_window = linspace(-4,4,size(posx_window,1))';
        speed_high_to_low_transition_pre(z,:) = speed2D(posx_window(1:size(posx_window,1)/2,1),posy_window(1:size(posx_window,1)/2,1),post_window(1:size(posx_window,1)/2,1));
        speed_high_to_low_transition_post(z,:) = speed2D(posx_window(size(posx_window,1)/2:end,1),posy_window(size(posx_window,1)/2:end,1),post_window(size(posx_window,1)/2:end,1));
    end
end

for z = 1:size(speed_high_to_low_transition_pre,1)
    if sum(speed_high_to_low_transition_pre(z,:)) == 0
        preremove(z) = false;
    else
        preremove(z) = true;
    end
end

for z = 1:size(speed_high_to_low_transition_post,1)
    if sum(speed_high_to_low_transition_post(z,:)) == 0
        postremove(z) = false;
    else
        postremove(z) = true;
    end
end

try
results.speed_high_to_low_transition_post = speed_high_to_low_transition_post(postremove,:);
results.speed_high_to_low_transition_pre = speed_high_to_low_transition_pre(preremove,:);

results.high_low_pretransition_mean = mean(results.speed_high_to_low_transition_pre,2);
results.high_low_postransition_mean = mean(results.speed_high_to_low_transition_post,2);
catch 
    keyboard
end

clear posx_window posy_window post_window preremove postremove
%%% CONTROL
pos_index_high_to_low_control = strfind(searchlog_control',[1 0]);

for z = 1:size(pos_index_high_to_low_control,2)
    if pos_index_high_to_low_control(1,z) - (50*4) <= 0 || pos_index_high_to_low_control(1,z) + (50*4) > size(xpos,1) %%% check to make sure they're no out of range
        continue
    else
        posx_window = xpos(pos_index_high_to_low_control(1,z)-(fs*4):pos_index_high_to_low_control(1,z)+(fs*4),1);
        posy_window = ypos(pos_index_high_to_low_control(1,z)-(fs*4):pos_index_high_to_low_control(1,z)+(fs*4),1);
        post_window = linspace(-4,4,size(posx_window,1))';
        speed_high_to_low_transition_pre_control(z,:) = speed2D(posx_window(1:size(posx_window,1)/2,1),posy_window(1:size(posx_window,1)/2,1),post_window(1:size(posx_window,1)/2,1));
        speed_high_to_low_transition_post_control(z,:) = speed2D(posx_window(size(posx_window,1)/2:end,1),posy_window(size(posx_window,1)/2:end,1),post_window(size(posx_window,1)/2:end,1));
    end
end

for z = 1:size(speed_high_to_low_transition_pre_control,1)
    if sum(speed_high_to_low_transition_pre_control(z,:)) == 0
        preremove(z) = false;
    else
        preremove(z) = true;
    end
end

for z = 1:size(speed_high_to_low_transition_post_control,1)
    if sum(speed_high_to_low_transition_post_control(z,:)) == 0
        postremove(z) = false;
    else
        postremove(z) = true;
    end
end

results.speed_high_to_low_transition_post_control = speed_high_to_low_transition_post_control(postremove,:);
results.speed_high_to_low_transition_pre_control = speed_high_to_low_transition_pre_control(preremove,:);

results.high_low_pretransition_mean_control = mean(results.speed_high_to_low_transition_pre_control,2);
results.high_low_postransition_mean_control = mean(results.speed_high_to_low_transition_post_control,2);

clear posx_window posy_window post_window preremove postremove