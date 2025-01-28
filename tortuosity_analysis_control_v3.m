function [results] = tortuosity_analysis_control_v3(tableaux,high_low_marker)
if high_low_marker == 0
%%% LOW TO HIGH TORTUOSITY TRANSITIONS
 tab = (tableaux);
 if isstruct(tab)
 tab = struct2table(tab);
 else
 end
%%% coherence at transition

 for i = 1:size(tab.coherence_wavelet_transition_low_high,1)
     for ii = 1:size(tab.coherence_wavelet_transition_low_high{i,1},2)     
       low_high_coherence(:,:,ii) = tab.coherence_wavelet_transition_low_high{i,1}{1,ii};
     end
     low_high.coherence_transition_low_high(:,:,i) = nanmean(low_high_coherence,3);
 end

 results.low_high.coherence_transition_low_high_by_animal = low_high.coherence_transition_low_high;
 results.low_high.coherence_transition_low_high = nanmean(low_high.coherence_transition_low_high,3);

%%% power at transition

for z = 1:size(tableaux,1)
    for zz = 1:size(tableaux.ch_one_transition_low_high_power_db{z,1},2)
        if isinf(sum(sum(tableaux.ch_one_transition_low_high_power_db{z,1}{zz})))
            day_ch_one(1:70,1:2001,zz) = NaN;
        else
        day_ch_one(:,:,zz) = tableaux.ch_one_transition_low_high_power_db{z,1}{zz};
        end

         if isinf(sum(sum(tableaux.ch_two_transition_low_high_power_db{z,1}{zz})))
            day_ch_two(1:70,1:2001,zz) = NaN;
        else
        day_ch_two(:,:,zz) = tableaux.ch_two_transition_low_high_power_db{z,1}{zz};
        end
        
    end
    channel_one(z,:,:) = squeeze(nanmean(day_ch_one,3));
    channel_two(z,:,:) = squeeze(nanmean(day_ch_two,3));
end

results.low_high.ch_one_power_transition_low_high = squeeze(nanmean(channel_one,1));
results.low_high.ch_two_power_transition_low_high = squeeze(nanmean(channel_two,1));

results.low_high.ch_one_power_by_animal = channel_one;
results.low_high.ch_one_power_by_animal = channel_two;

clear channel_one channel_two day_ch_one day_ch_two
%%% continuous angle difference for 1000ms pre and 1000ms post change
%%% (hilbert determined theta-filtered angle)

for z = 1:size(tab.angle_difference_pre,1)
    %%% average over transitions in each recording
       for zz = 1:size(tab.angle_difference_pre{z,1},2)
        if z == 1
             results.low_high.angle_prechange_session_low_high(1:1000,z) = tab.angle_difference_pre{z,1}{1,zz};
             results.low_high.angle_postchange_session_low_high(1:1001,z) = tab.angle_difference_post{z,1}{1,zz};
        else
        
             results.low_high.angle_prechange_session_low_high(1:1000,z+zz) = tab.angle_difference_pre{z,1}{1,zz};
             results.low_high.angle_postchange_session_low_high(1:1001,z+zz) = tab.angle_difference_post{z,1}{1,zz};
          
        end
       end
end

for e = 1:size(results.low_high.angle_prechange_session_low_high,2)
    if sum(results.low_high.angle_prechange_session_low_high(:,e)) == 0
angle_delete_pre(:,e) = 0;
    else
        angle_delete_pre(:,e) = 1;
    end
end

for e = 1:size(results.low_high.angle_postchange_session_low_high,2)
    if sum(results.low_high.angle_postchange_session_low_high(:,e)) == 0
angle_delete_post(:,e) = 0;
    else
        angle_delete_post(:,e) = 1;
    end
end

results.low_high.angle_prechange_session_low_high(:,~angle_delete_pre) = [];
results.low_high.angle_postchange_session_low_high(:,~angle_delete_post) = [];

for e = 1:size(results.low_high.angle_prechange_session_low_high,1)
results.low_high.mean_angle_prechange_session_low_high(e,:)= circ_mean(results.low_high.angle_prechange_session_low_high(e,:));
end
for e = 1:size(results.low_high.angle_postchange_session_low_high,1)
results.low_high.mean_angle_postchange_session_low_high(e,:) = circ_mean(results.low_high.angle_postchange_session_low_high(e,:));
end

%%% peak offsets (ms)
%%% channel one to channel two

one_two_prechange_peak_offset_seconds = cat(2,tab.one_two_prechange_peak_offset_seconds);
one_two_postchange_peak_offset_seconds = cat(2,tab.one_two_postchange_peak_offset_seconds);

for z = 1:size(one_two_prechange_peak_offset_seconds,1)
  results.low_high.one_two_prechange_mean_peak_offset_in_seconds(z,:) = nanmean(cell2mat(one_two_prechange_peak_offset_seconds{z}));
  results.low_high.one_two_prechange_sterr_peak_offset_in_seconds(z,:) = nanstd(cell2mat(one_two_prechange_peak_offset_seconds{z}))/sqrt(size(one_two_prechange_peak_offset_seconds{z},2));

  results.low_high.one_two_postchange_mean_peak_offset_in_seconds(z,:) = nanmean(cell2mat(one_two_postchange_peak_offset_seconds{z}));
  results.low_high.one_two_postchange_sterr_peak_offset_in_seconds(z,:) = nanstd(cell2mat(one_two_postchange_peak_offset_seconds{z}))/sqrt(size(one_two_postchange_peak_offset_seconds{z},2));
end
    
one_two_prechange_peak_offset_angle = cat(2,tab.one_two_prechange_peak_angle_offset);
one_two_postchange_peak_offset_angle = cat(2,tab.one_two_postchange_peak_angle_offset);

for z = 1:size(one_two_prechange_peak_offset_seconds,1)
  results.low_high.one_two_mean_peak_offset_angle(z,:) = circ_mean(rmmissing(cell2mat(one_two_prechange_peak_offset_angle{z}')));
  results.low_high.one_two_sterr_peak_offset_angle(z,:) = circ_std(rmmissing(cell2mat(one_two_prechange_peak_offset_angle{z}')))/sqrt(size(one_two_prechange_peak_offset_angle{z},2));
end

%%% peak offsets (ms)
%%% channel two to channel one

two_one_prechange_peak_offset_seconds = cat(2,tab.two_one_prechange_peak_offset_seconds);
two_one_postchange_peak_offset_seconds = cat(2,tab.two_one_postchange_peak_offset_seconds);

for z = 1:size(two_one_prechange_peak_offset_seconds,1)
  results.low_high.two_one_prechange_mean_peak_offset_in_seconds(z,:) = nanmean(cell2mat(two_one_prechange_peak_offset_seconds{z}));
  results.low_high.two_one_prechange_sterr_peak_offset_in_seconds(z,:) = nanstd(cell2mat(two_one_prechange_peak_offset_seconds{z}))/sqrt(size(two_one_prechange_peak_offset_seconds{z},2));

  results.low_high.two_one_postchange_mean_peak_offset_in_seconds(z,:) = nanmean(cell2mat(two_one_postchange_peak_offset_seconds{z}));
  results.low_high.two_one_postchange_sterr_peak_offset_in_seconds(z,:) = nanstd(cell2mat(two_one_postchange_peak_offset_seconds{z}))/sqrt(size(two_one_postchange_peak_offset_seconds{z},2));
end
    
two_one_prechange_peak_offset_angle = cat(2,tab.two_one_prechange_peak_angle_offset);
two_one_postchange_peak_offset_angle = cat(2,tab.two_one_postchange_peak_angle_offset);

for z = 1:size(two_one_prechange_peak_offset_seconds,1)
  results.low_high.two_one_mean_peak_offset_angle(z,:) = circ_mean(rmmissing(cell2mat(two_one_prechange_peak_offset_angle{z}')));
  results.low_high.two_one_sterr_peak_offset_angle(z,:) = circ_std(rmmissing(cell2mat(two_one_prechange_peak_offset_angle{z}')))/sqrt(size(two_one_prechange_peak_offset_angle{z},2));
end

elseif high_low_marker == 1
%% HIGH TO LOW TORTUOSITY TRANSITIONS
tab = (tableaux);
if isstruct(tab)
tab = struct2table(tab);
else
end
 results.high_low.coherence_transition_high_low = zeros(84,2001);
 for i = 1:size(tab.coherence_wavelet_transition_high_low,1)
     for ii = 1:size(tab.coherence_wavelet_transition_high_low{i,1},2)     
       high_low_coherence(:,:,ii) = tab.coherence_wavelet_transition_high_low{i,1}{1,ii};
     end
     high_low.coherence_transition_high_low(:,:,i) = nanmean(high_low_coherence,3);
 end

 results.high_low.coherence_transition_high_low_by_animal = high_low.coherence_transition_high_low;
 results.high_low.coherence_transition_high_low = nanmean(high_low.coherence_transition_high_low,3);

%%% power at transition
for z = 1:size(tableaux,1)
    for zz = 1:size(tableaux.ch_one_transition_high_low_power_db{z,1},2)
        if isinf(sum(sum(tableaux.ch_one_transition_high_low_power_db{z,1}{zz})))
            day_ch_one(1:70,1:2001,zz) = NaN;
        else
        day_ch_one(:,:,zz) = tableaux.ch_one_transition_high_low_power_db{z,1}{zz};
        end

         if isinf(sum(sum(tableaux.ch_two_transition_high_low_power_db{z,1}{zz})))
            day_ch_two(1:70,1:2001,zz) = NaN;
        else
        day_ch_two(:,:,zz) = tableaux.ch_two_transition_high_low_power_db{z,1}{zz};
        end
        
    end
    channel_one(z,:,:) = squeeze(nanmean(day_ch_one,3));
    channel_two(z,:,:) = squeeze(nanmean(day_ch_two,3));
end

results.high_low.ch_one_power_transition_high_low = squeeze(nanmean(channel_one,1));
results.high_low.ch_two_power_transition_high_low = squeeze(nanmean(channel_two,1));

results.high_low.ch_one_power_by_animal = channel_one;
results.high_low.ch_two_power_by_animal = channel_two;

clear channel_one channel_two day_ch_one day_ch_two
%%% continuous angle difference for 1000ms pre and 1000ms post change
%%% (hilbert determined theta-filtered angle)

for z = 1:size(tab.angle_difference_pre,1)
    %%% average over transitions in each recording
       for zz = 1:size(tab.angle_difference_pre{z,1},2)
        if z == 1
             results.high_low.angle_prechange_session_high_low(1:1000,z) = tab.angle_difference_pre{z,1}{1,zz};
             results.high_low.angle_postchange_session_high_low(1:1001,z) = tab.angle_difference_post{z,1}{1,zz};
        else
        
             results.high_low.angle_prechange_session_high_low(1:1000,z+zz) = tab.angle_difference_pre{z,1}{1,zz};
             results.high_low.angle_postchange_session_high_low(1:1001,z+zz) = tab.angle_difference_post{z,1}{1,zz};
          
        end
       end
end

for e = 1:size(results.high_low.angle_prechange_session_high_low,2)
    if sum(results.high_low.angle_prechange_session_high_low(:,e)) == 0
angle_delete_pre(:,e) = 0;
    else
        angle_delete_pre(:,e) = 1;
    end
end

for e = 1:size(results.high_low.angle_postchange_session_high_low,2)
    if sum(results.high_low.angle_postchange_session_high_low(:,e)) == 0
angle_delete_post(:,e) = 0;
    else
        angle_delete_post(:,e) = 1;
    end
end

results.high_low.angle_prechange_session_high_low(:,~angle_delete_pre) = [];
results.high_low.angle_postchange_session_high_low(:,~angle_delete_post) = [];

for e = 1:size(results.high_low.angle_prechange_session_high_low,1)
results.high_low.mean_angle_prechange_session_high_low(e,:)= circ_mean(results.high_low.angle_prechange_session_high_low(e,:));
end
for e = 1:size(results.high_low.angle_postchange_session_high_low,1)
results.high_low.mean_angle_postchange_session_high_low(e,:) = circ_mean(results.high_low.angle_postchange_session_high_low(e,:));
end


%%% peak offsets (ms)
%%% channel one to channel two

one_two_prechange_peak_offset_seconds = cat(2,tab.one_two_prechange_peak_offset_seconds);
one_two_postchange_peak_offset_seconds = cat(2,tab.one_two_postchange_peak_offset_seconds);

for z = 1:size(one_two_prechange_peak_offset_seconds,1)
  results.high_low.one_two_prechange_mean_peak_offset_in_seconds(z,:) = nanmean(cell2mat(one_two_prechange_peak_offset_seconds{z}));
  results.high_low.one_two_prechange_sterr_peak_offset_in_seconds(z,:) = nanstd(cell2mat(one_two_prechange_peak_offset_seconds{z}))/sqrt(size(one_two_prechange_peak_offset_seconds{z},2));

  results.high_low.one_two_postchange_mean_peak_offset_in_seconds(z,:) = nanmean(cell2mat(one_two_postchange_peak_offset_seconds{z}));
  results.high_low.one_two_postchange_sterr_peak_offset_in_seconds(z,:) = nanstd(cell2mat(one_two_postchange_peak_offset_seconds{z}))/sqrt(size(one_two_postchange_peak_offset_seconds{z},2));
end
    
one_two_prechange_peak_offset_angle = cat(2,tab.one_two_prechange_peak_angle_offset);
one_two_postchange_peak_offset_angle = cat(2,tab.one_two_postchange_peak_angle_offset);

for z = 1:size(one_two_prechange_peak_offset_seconds,1)
  results.high_low.one_two_mean_peak_offset_angle(z,:) = circ_mean(rmmissing(cell2mat(one_two_prechange_peak_offset_angle{z}')));
  results.high_low.one_two_sterr_peak_offset_angle(z,:) = circ_std(rmmissing(cell2mat(one_two_prechange_peak_offset_angle{z}')))/sqrt(size(one_two_prechange_peak_offset_angle{z},2));
end

%%% peak offsets (ms)
%%% channel two to channel one

two_one_prechange_peak_offset_seconds = cat(2,tab.two_one_prechange_peak_offset_seconds);
two_one_postchange_peak_offset_seconds = cat(2,tab.two_one_postchange_peak_offset_seconds);

for z = 1:size(two_one_prechange_peak_offset_seconds,1)
  results.high_low.two_one_prechange_mean_peak_offset_in_seconds(z,:) = nanmean(cell2mat(two_one_prechange_peak_offset_seconds{z}));
  results.high_low.two_one_prechange_sterr_peak_offset_in_seconds(z,:) = nanstd(cell2mat(two_one_prechange_peak_offset_seconds{z}))/sqrt(size(two_one_prechange_peak_offset_seconds{z},2));

  results.high_low.two_one_postchange_mean_peak_offset_in_seconds(z,:) = nanmean(cell2mat(two_one_postchange_peak_offset_seconds{z}));
  results.high_low.two_one_postchange_sterr_peak_offset_in_seconds(z,:) = nanstd(cell2mat(two_one_postchange_peak_offset_seconds{z}))/sqrt(size(two_one_postchange_peak_offset_seconds{z},2));
end
    
two_one_prechange_peak_offset_angle = cat(2,tab.two_one_prechange_peak_angle_offset);
two_one_postchange_peak_offset_angle = cat(2,tab.two_one_postchange_peak_angle_offset);

for z = 1:size(two_one_prechange_peak_offset_seconds,1)
  results.high_low.two_one_mean_peak_offset_angle(z,:) = circ_mean(rmmissing(cell2mat(two_one_prechange_peak_offset_angle{z}')));
  results.high_low.two_one_sterr_peak_offset_angle(z,:) = circ_std(rmmissing(cell2mat(two_one_prechange_peak_offset_angle{z}')))/sqrt(size(two_one_prechange_peak_offset_angle{z},2));
end

end