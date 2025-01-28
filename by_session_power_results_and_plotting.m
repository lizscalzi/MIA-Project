function [freq_power_stats] = by_session_power_results_and_plotting(ctl,ctl_ctl_segs,mia,mia_ctl_segs,plotfig,freq,time)
keyboard
theta_freqs = freq(6:11);
%% CHANNEL ONE
%% low to high
%%% around transitions

ctl_mean_around_transitions_low_high = mean(ctl.around_transitions.low_high_peri_channel_one_freq_domain,2);
ctl_sem_around_transitions_low_high = std(ctl.around_transitions.low_high_peri_channel_one_freq_domain,0,2)./sqrt(size(ctl.around_transitions.low_high_peri_channel_one_freq_domain,2));
mia_mean_around_transitions_low_high = mean(mia.around_transitions.low_high_peri_channel_one_freq_domain,2);
mia_sem_around_transitions_low_high = std(mia.around_transitions.low_high_peri_channel_one_freq_domain,0,2)./sqrt(size(ctl.around_transitions.low_high_peri_channel_one_freq_domain,2));
ctl_ctl_mean_around_transitions_low_high = mean(ctl_ctl_segs.around_transitions.low_high_peri_channel_one_freq_domain,2);
ctl_ctl_sem_around_transitions_low_high = std(ctl_ctl_segs.around_transitions.low_high_peri_channel_one_freq_domain,0,2)./sqrt(size(ctl.around_transitions.low_high_peri_channel_one_freq_domain,2));
mia_ctl_mean_around_transitions_low_high = mean(mia_ctl_segs.around_transitions.low_high_peri_channel_one_freq_domain,2);
mia_ctl_sem_around_transitions_low_high = std(mia_ctl_segs.around_transitions.low_high_peri_channel_one_freq_domain,0,2)./sqrt(size(ctl.around_transitions.low_high_peri_channel_one_freq_domain,2));

if plotfig == true
    figure(1)
    shadedErrorBar(theta_freqs,ctl_mean_around_transitions_low_high,ctl_sem_around_transitions_low_high,'lineProps','k-');
    hold on
    shadedErrorBar(theta_freqs,ctl_ctl_mean_around_transitions_low_high,ctl_ctl_sem_around_transitions_low_high,'lineProps','r-');
    title('Theta power during low to high transisitions, control animals');
    ylabel('Power difference from baseline (dB)');
    xlabel('Frequency (Hz)');

    figure(2)
    shadedErrorBar(theta_freqs,mia_mean_around_transitions_low_high,ctl_sem_around_transitions_low_high,'lineProps','k-');
    hold on
    shadedErrorBar(theta_freqs,mia_ctl_mean_around_transitions_low_high,ctl_ctl_sem_around_transitions_low_high,'lineProps','r-');
    title('Theta power during low to high transisitions, MIA animals');
    ylabel('Power difference from baseline (dB)');
    xlabel('Frequency (Hz)');
end

%% high to low
ctl_mean_around_transitions_high_low = mean(ctl.around_transitions.high_low_peri_channel_one_freq_domain,2);
ctl_sem_around_transitions_high_low = std(ctl.around_transitions.high_low_peri_channel_one_freq_domain,0,2)./sqrt(size(ctl.around_transitions.high_low_peri_channel_one_freq_domain,2));
mia_mean_around_transitions_high_low = mean(mia.around_transitions.high_low_peri_channel_one_freq_domain,2);
mia_sem_around_transitions_high_low = std(mia.around_transitions.high_low_peri_channel_one_freq_domain,0,2)./sqrt(size(ctl.around_transitions.high_low_peri_channel_one_freq_domain,2));
ctl_ctl_mean_around_transitions_high_low = mean(ctl_ctl_segs.around_transitions.high_low_peri_channel_one_freq_domain,2);
ctl_ctl_sem_around_transitions_high_low = std(ctl_ctl_segs.around_transitions.high_low_peri_channel_one_freq_domain,0,2)./sqrt(size(ctl.around_transitions.high_low_peri_channel_one_freq_domain,2));
mia_ctl_mean_around_transitions_high_low = mean(mia_ctl_segs.around_transitions.high_low_peri_channel_one_freq_domain,2);
mia_ctl_sem_around_transitions_high_low = std(mia_ctl_segs.around_transitions.high_low_peri_channel_one_freq_domain,0,2)./sqrt(size(ctl.around_transitions.high_low_peri_channel_one_freq_domain,2));

if plotfig == true
    figure(2)
    shadedErrorBar(theta_freqs,ctl_mean_around_transitions_high_low,ctl_sem_around_transitions_high_low,'lineProps','k-');
    hold on
    shadedErrorBar(theta_freqs,ctl_ctl_mean_around_transitions_high_low,ctl_ctl_sem_around_transitions_high_low,'lineProps','r-');
    title('Theta power during high to low transisitions, control animals');
    ylabel('Power difference from baseline (dB)');
    xlabel('Frequency (Hz)');

    figure(2)
    shadedErrorBar(theta_freqs,mia_mean_around_transitions_high_low,ctl_sem_around_transitions_high_low,'lineProps','k-');
    hold on
    shadedErrorBar(theta_freqs,mia_ctl_mean_around_transitions_high_low,ctl_ctl_sem_around_transitions_high_low,'lineProps','r-');
    title('Theta power during high to low transisitions, MIA animals');
    ylabel('Power difference from baseline (dB)');
    xlabel('Frequency (Hz)');
end



%%%% time domain
%%% low to high

around_transitions_time = linspace(-1,1,251);
ctl_mean_around_transitions_low_high_time = mean(ctl.around_transitions.low_high_peri_channel_one_time,1);
ctl_sem_around_transitions_low_high_time = std(ctl.around_transitions.low_high_peri_channel_one_time,0,1)./sqrt(size(ctl.around_transitions.low_high_peri_channel_one_time,1));

mia_mean_around_transitions_low_high_time = mean(mia.around_transitions.low_high_peri_channel_one_time,1);
mia_sem_around_transitions_low_high_time = std(mia.around_transitions.low_high_peri_channel_one_time,0,1)./sqrt(size(ctl.around_transitions.low_high_peri_channel_one_time,1));


ctl_ctl_mean_around_transitions_low_high_time = mean(ctl_ctl_segs.around_transitions.low_high_peri_channel_one_time,1);
ctl_ctl_sem_around_transitions_low_high_time = std(ctl_ctl_segs.around_transitions.low_high_peri_channel_one_time,0,1)./sqrt(size(ctl.around_transitions.low_high_peri_channel_one_time,1));

mia_ctl_mean_around_transitions_low_high_time = mean(mia_ctl_segs.around_transitions.low_high_peri_channel_one_time,1);
mia_ctl_sem_around_transitions_low_high_time = std(mia_ctl_segs.around_transitions.low_high_peri_channel_one_time,0,1)./sqrt(size(ctl.around_transitions.low_high_peri_channel_one_time,1));

if plotfig == true
    figure(1)
    subplot(1,2,1)
    shadedErrorBar(around_transitions_time,ctl_mean_around_transitions_low_high_time,ctl_sem_around_transitions_low_high_time,'lineProps','k-');
    hold on
    shadedErrorBar(around_transitions_time,ctl_ctl_mean_around_transitions_low_high_time,ctl_ctl_sem_around_transitions_low_high_time,'lineProps','r-');
    title('Spectral power, control animals');
    subtitle('Peri Transition, Low to High')
    ylabel('Power difference from baseline (dB)');
    legend('Transition segments','  ',' ',' ','Control Segments',' ');
    legend('boxoff');
    xline(0,'k--');
    xlabel('Time (s)');
    subplot(1,2,2)
    shadedErrorBar(around_transitions_time',mia_mean_around_transitions_low_high_time,ctl_sem_around_transitions_low_high_time,'lineProps','k-');
    hold on
    shadedErrorBar(around_transitions_time',mia_ctl_mean_around_transitions_low_high_time,ctl_ctl_sem_around_transitions_low_high_time,'lineProps','r-');
    title('Spectral power, MIA animals');
    subtitle('Peri Transition, Low to High')
    ylabel('Power difference from baseline (dB)');
    xlabel('Time (s)');
    legend('Transition segments','  ',' ',' ','Control Segments',' ');
    legend('boxoff');
     xline(0,'k--');
end




%% CHANNEL TWO
%% low to high

%%% around transitions

ctl_mean_around_transitions_low_high = mean(ctl.around_transitions.low_high_peri_channel_two_freq_domain,2);
ctl_sem_around_transitions_low_high = std(ctl.around_transitions.low_high_peri_channel_two_freq_domain,0,2)./sqrt(size(ctl.around_transitions.low_high_peri_channel_two_freq_domain,2));

mia_mean_around_transitions_low_high = mean(mia.around_transitions.low_high_peri_channel_two_freq_domain,2);
mia_sem_around_transitions_low_high = std(mia.around_transitions.low_high_peri_channel_two_freq_domain,0,2)./sqrt(size(ctl.around_transitions.low_high_peri_channel_two_freq_domain,2));


ctl_ctl_mean_around_transitions_low_high = mean(ctl_ctl_segs.around_transitions.low_high_peri_channel_two_freq_domain,2);
ctl_ctl_sem_around_transitions_low_high = std(ctl_ctl_segs.around_transitions.low_high_peri_channel_two_freq_domain,0,2)./sqrt(size(ctl.around_transitions.low_high_peri_channel_two_freq_domain,2));

mia_ctl_mean_around_transitions_low_high = mean(mia_ctl_segs.around_transitions.low_high_peri_channel_two_freq_domain,2);
mia_ctl_sem_around_transitions_low_high = std(mia_ctl_segs.around_transitions.low_high_peri_channel_two_freq_domain,0,2)./sqrt(size(ctl.around_transitions.low_high_peri_channel_two_freq_domain,2));

if plotfig == true
    figure(1)
    shadedErrorBar(theta_freqs,ctl_mean_around_transitions_low_high,ctl_sem_around_transitions_low_high,'lineProps','k-');
    hold on
    shadedErrorBar(theta_freqs,ctl_ctl_mean_around_transitions_low_high,ctl_ctl_sem_around_transitions_low_high,'lineProps','r-');
    title('Theta power during low to high transisitions, control animals');
    ylabel('Power difference from baseline (dB)');
    xlabel('Frequency (Hz)');

    figure(2)
    shadedErrorBar(theta_freqs,mia_mean_around_transitions_low_high,ctl_sem_around_transitions_low_high,'lineProps','k-');
    hold on
    shadedErrorBar(theta_freqs,mia_ctl_mean_around_transitions_low_high,ctl_ctl_sem_around_transitions_low_high,'lineProps','r-');
    title('Theta power during low to high transisitions, MIA animals');
    ylabel('Power difference from baseline (dB)');
    xlabel('Frequency (Hz)');
end

%% high to low
ctl_mean_around_transitions_high_low = mean(ctl.around_transitions.high_low_peri_channel_two_freq_domain,2);
ctl_sem_around_transitions_high_low = std(ctl.around_transitions.high_low_peri_channel_two_freq_domain,0,2)./sqrt(size(ctl.around_transitions.high_low_peri_channel_two_freq_domain,2));

mia_mean_around_transitions_high_low = mean(mia.around_transitions.high_low_peri_channel_two_freq_domain,2);
mia_sem_around_transitions_high_low = std(mia.around_transitions.high_low_peri_channel_two_freq_domain,0,2)./sqrt(size(ctl.around_transitions.high_low_peri_channel_two_freq_domain,2));


ctl_ctl_mean_around_transitions_high_low = mean(ctl_ctl_segs.around_transitions.high_low_peri_channel_two_freq_domain,2);
ctl_ctl_sem_around_transitions_high_low = std(ctl_ctl_segs.around_transitions.high_low_peri_channel_two_freq_domain,0,2)./sqrt(size(ctl.around_transitions.high_low_peri_channel_two_freq_domain,2));

mia_ctl_mean_around_transitions_high_low = mean(mia_ctl_segs.around_transitions.high_low_peri_channel_two_freq_domain,2);
mia_ctl_sem_around_transitions_high_low = std(mia_ctl_segs.around_transitions.high_low_peri_channel_two_freq_domain,0,2)./sqrt(size(ctl.around_transitions.high_low_peri_channel_two_freq_domain,2));

if plotfig == true
    figure(1)
    shadedErrorBar(theta_freqs,ctl_mean_around_transitions_high_low,ctl_sem_around_transitions_high_low,'lineProps','k-');
    hold on
    shadedErrorBar(theta_freqs,ctl_ctl_mean_around_transitions_high_low,ctl_ctl_sem_around_transitions_high_low,'lineProps','r-');
    title('Theta power during high to low transisitions, control animals');
    ylabel('Power difference from baseline (dB)');
    xlabel('Frequency (Hz)');

    figure(2)
    shadedErrorBar(theta_freqs,mia_mean_around_transitions_high_low,ctl_sem_around_transitions_high_low,'lineProps','k-');
    hold on
    shadedErrorBar(theta_freqs,mia_ctl_mean_around_transitions_high_low,ctl_ctl_sem_around_transitions_high_low,'lineProps','r-');
    title('Theta power during high to low transisitions, MIA animals');
    ylabel('Power difference from baseline (dB)');
    xlabel('Frequency (Hz)');
end