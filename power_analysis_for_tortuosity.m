function [results] = power_analysis_for_tortuosity(lh_one,lh_two,hl_one,hl_two,frequencies,times)

index_to_one_second_around_transitions = 876:1126;
index_to_two_to_one_before_transitions = 501:751;
index_to_one_to_two_post_transitions = 1251:1501;
index_to_theta = 6:11;



%% around tranisitions
for i = 1:size(lh_one,1)
    lh_peri_change_power_ch_one(i,:,:) = squeeze(lh_one(i,index_to_theta,index_to_one_second_around_transitions));
    lh_peri_change_power_ch_two(i,:,:) = squeeze(lh_two(i,index_to_theta,index_to_one_second_around_transitions));
    hl_peri_change_power_ch_one(i,:,:) = squeeze(hl_one(i,index_to_theta,index_to_one_second_around_transitions));
    hl_peri_change_power_ch_two(i,:,:) = squeeze(hl_two(i,index_to_theta,index_to_one_second_around_transitions));
end

for i = 1:size(lh_peri_change_power_ch_one,1)
    results.around_transitions.low_high_peri_channel_one_time(i,:) = mean(lh_peri_change_power_ch_one(i,:,:),2);
    results.around_transitions.low_high_peri_channel_two_time(i,:) = mean(lh_peri_change_power_ch_two(i,:,:),2);
    results.around_transitions.high_low_peri_channel_one_time(i,:) =mean(hl_peri_change_power_ch_one(i,:,:),2);
    results.around_transitions.high_low_peri_channel_two_time(i,:) = mean(hl_peri_change_power_ch_two(i,:,:),2);
end


for i = 1:size(lh_peri_change_power_ch_one,1)
    results.around_transitions.low_high_peri_channel_one_freq_domain(:,i) = mean(lh_peri_change_power_ch_one(i,:,:),3);
    results.around_transitions.low_high_peri_channel_two_freq_domain(:,i) = mean(lh_peri_change_power_ch_two(i,:,:),3);
    results.around_transitions.high_low_peri_channel_one_freq_domain(:,i) =mean(hl_peri_change_power_ch_one(i,:,:),3);
    results.around_transitions.high_low_peri_channel_two_freq_domain(:,i) = mean(hl_peri_change_power_ch_two(i,:,:),3);
end
clear lh_peri_change_power_ch_one lh_peri_change_power_ch_two hl_peri_change_power_ch_one hl_peri_change_power_ch_two
%% -2 to 1 second prior to transition
for i = 1:size(lh_one,1)
    lh_peri_change_power_ch_one(i,:,:) = squeeze(lh_one(i,index_to_theta,index_to_two_to_one_before_transitions));
    lh_peri_change_power_ch_two(i,:,:) = squeeze(lh_two(i,index_to_theta,index_to_two_to_one_before_transitions));
    hl_peri_change_power_ch_one(i,:,:) = squeeze(hl_one(i,index_to_theta,index_to_two_to_one_before_transitions));
    hl_peri_change_power_ch_two(i,:,:) = squeeze(hl_two(i,index_to_theta,index_to_two_to_one_before_transitions));
end

for i = 1:size(lh_peri_change_power_ch_one,1)
    results.minus_two_to_one_prior.low_high_peri_channel_one_time(i,:) = mean(lh_peri_change_power_ch_one(i,:,:),2);
    results.minus_two_to_one_prior.low_high_peri_channel_two_time(i,:) = mean(lh_peri_change_power_ch_two(i,:,:),2);
    results.minus_two_to_one_prior.high_low_peri_channel_one_time(i,:) =mean(hl_peri_change_power_ch_one(i,:,:),2);
    results.minus_two_to_one_prior.high_low_peri_channel_two_time(i,:) = mean(hl_peri_change_power_ch_two(i,:,:),2);
end


for i = 1:size(lh_peri_change_power_ch_one,1)
    results.minus_two_to_one_prior.low_high_peri_channel_one_freq_domain(:,i) = mean(lh_peri_change_power_ch_one(i,:,:),3);
    results.minus_two_to_one_prior.low_high_peri_channel_two_freq_domain(:,i) = mean(lh_peri_change_power_ch_two(i,:,:),3);
    results.minus_two_to_one_prior.high_low_peri_channel_one_freq_domain(:,i) =mean(hl_peri_change_power_ch_one(i,:,:),3);
    results.minus_two_to_one_prior.high_low_peri_channel_two_freq_domain(:,i) = mean(hl_peri_change_power_ch_two(i,:,:),3);
end
clear lh_peri_change_power_ch_one lh_peri_change_power_ch_two hl_peri_change_power_ch_one hl_peri_change_power_ch_two

%% 1 to 2 seconds post transition
for i = 1:size(lh_one,1)
    lh_peri_change_power_ch_one(i,:,:) = squeeze(lh_one(i,index_to_theta,index_to_one_to_two_post_transitions));
    lh_peri_change_power_ch_two(i,:,:) = squeeze(lh_two(i,index_to_theta,index_to_one_to_two_post_transitions));
    hl_peri_change_power_ch_one(i,:,:) = squeeze(hl_one(i,index_to_theta,index_to_one_to_two_post_transitions));
    hl_peri_change_power_ch_two(i,:,:) = squeeze(hl_two(i,index_to_theta,index_to_one_to_two_post_transitions));
end

for i = 1:size(lh_peri_change_power_ch_one,1)
    results.one_to_two_seconds_post.low_high_peri_channel_one_time(i,:) = mean(lh_peri_change_power_ch_one(i,:,:),2);
    results.one_to_two_seconds_post.low_high_peri_channel_two_time(i,:) = mean(lh_peri_change_power_ch_two(i,:,:),2);
    results.one_to_two_seconds_post.high_low_peri_channel_one_time(i,:) =mean(hl_peri_change_power_ch_one(i,:,:),2);
    results.one_to_two_seconds_post.high_low_peri_channel_two_time(i,:) = mean(hl_peri_change_power_ch_two(i,:,:),2);
end


for i = 1:size(lh_peri_change_power_ch_one,1)
    results.minus_two_to_one_prior.low_high_peri_channel_one_freq_domain(:,i) = mean(lh_peri_change_power_ch_one(i,:,:),3);
    results.minus_two_to_one_prior.low_high_peri_channel_two_freq_domain(:,i) = mean(lh_peri_change_power_ch_two(i,:,:),3);
    results.minus_two_to_one_prior.high_low_peri_channel_one_freq_domain(:,i) =mean(hl_peri_change_power_ch_one(i,:,:),3);
    results.minus_two_to_one_prior.high_low_peri_channel_two_freq_domain(:,i) = mean(hl_peri_change_power_ch_two(i,:,:),3);
end

clear lh_peri_change_power_ch_one lh_peri_change_power_ch_two hl_peri_change_power_ch_one hl_peri_change_power_ch_two