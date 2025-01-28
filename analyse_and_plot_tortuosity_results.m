clear all

%% Do lfp tortuosity analysis
%% CONTROL
ctl_lfp_file = uigetfile();
load(ctl_lfp_file);
%load('/Volumes/1TB ALPHA/MIA_Coherence_Results/CTL_Tortuosity_LFP_Results_low_high_v3.mat');
tableaux_ctl_low_high = struct2table(tortuosity_low_high_ctl); clear tortuosity_low_high_ctl
[tortuosity_control_low_high] = tortuosity_analysis_control_v3(tableaux_ctl_low_high,0);
clear tableaux_ctl_low_high

%% MIA
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/MIA_Tortuosity_LFP_Results_low_high_v3.mat');
tableaux_mia_low_high = struct2table(tortuosity_low_high_MIA); clear tortuosity_low_high_mia
[tortuosity_mia_low_high] = tortuosity_analysis_control_v3(tableaux_mia_low_high,0);
clear tableaux_mia_low_high
%% CONTROL
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/CTL_Tortuosity_LFP_Results_high_low_v3.mat');
tableaux_ctl_high_low = struct2table(tortuosity_high_low_ctl); clear tortuosity_high_low_ctl
[tortuosity_control_high_low] = tortuosity_analysis_control_v3(tableaux_ctl_high_low,1);
clear tableaux_ctl_high_low

%% MIA
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/MIA_Tortuosity_LFP_Results_high_low_v3.mat');
tableaux_mia_high_low = struct2table(tortuosity_high_low_MIA); clear tortuosity_high_low_mia
[tortuosity_mia_high_low] = tortuosity_analysis_control_v3(tableaux_mia_high_low,1);
clear tableaux_mia_high_low


load('/Volumes/1TB ALPHA/MIA_Coherence_Results/tortuosity_frequencies.mat');
coherence_times = tortuosity_freqs.transition_timebase;
coherence_freqs = tortuosity_freqs.coherence_freqs;
power_freqs = tortuosity_freqs.power_freqs;
power_times = tortuosity_freqs.transition_timebase;



%% FIGURES
%%%% coherence figure
figure(1)
set(0, 'DefaultFigureRenderer', 'painters');
subplot(2,2,1)
imagesc(coherence_times,coherence_freqs,tortuosity_control_low_high.low_high.coherence_transition_low_high,[0.3 0.6]);
set(gca,'YDir','normal')
xlabel('Low           Time from tortuosity switch (seconds)         High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Coherence");
title('Control, low to high');
subplot(2,2,2);
imagesc(coherence_times,coherence_freqs,tortuosity_mia_low_high.low_high.coherence_transition_low_high,[0.3 0.6]);
set(gca,'YDir','normal');
xlabel('Low           Time from tortuosity switch (seconds)        High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Coherence");
title('MIA, low to high');
subplot(2,2,3)
imagesc(coherence_times,coherence_freqs,tortuosity_control_high_low.high_low.coherence_transition_high_low,[0.3 0.6]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Coherence");
title('Control, high to low');
subplot(2,2,4)
imagesc(coherence_times,coherence_freqs,tortuosity_mia_high_low.high_low.coherence_transition_high_low,[0.3 0.6]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Coherence");
title('MIA, high to low');
sgtitle('Spectral coherence at tortuosity transitions');
colormap jet;







%%%% power figure
figure(2)
set(0, 'DefaultFigureRenderer', 'painters');
subplot(2,2,1)
contourf(power_times,power_freqs,imgaussfilt(tortuosity_control_low_high.low_high.ch_one_power_transition_low_high,4),'linecolor','none');
caxis([10 28]);
xlabel('Low           Time from tortuosity switch (seconds)         High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Power (dB)");
title('Control channel 1 (PFC)');
subplot(2,2,2);
contourf(power_times,power_freqs,imgaussfilt(tortuosity_mia_low_high.low_high.ch_one_power_transition_low_high,4),'linecolor','none');;
caxis([10 28]);
xlabel('Low           Time from tortuosity switch (seconds)        High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Power (dB)");
title('MIA channel 1 (PFC)');
subplot(2,2,3)
contourf(power_times,power_freqs,imgaussfilt(tortuosity_control_low_high.low_high.ch_two_power_transition_low_high,4),'linecolor','none');
caxis([10 28]);
xlabel('Low           Time from tortuosity switch (seconds)         High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Power (dB)");
title('Control channel 2 (HPC)');
subplot(2,2,4);
contourf(power_times,power_freqs,imgaussfilt(tortuosity_mia_low_high.low_high.ch_two_power_transition_low_high,4),'linecolor','none');
caxis([10 28]);
xlabel('Low           Time from tortuosity switch (seconds)        High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Power (dB)");
title('MIA Channel 2 (HPC)');
sgtitle('Spectral power during low to high tortuosity transitions')
colormap jet;

figure(3)
set(0, 'DefaultFigureRenderer', 'painters');
subplot(2,2,1)
contourf(power_times,power_freqs,imgaussfilt(tortuosity_control_high_low.high_low.ch_one_power_transition_high_low,4),'linecolor','none');
caxis([10 28]);
xlabel('High          Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Power (dB)");
title('Control channel 1 (PFC)');
subplot(2,2,2);
contourf(power_times,power_freqs,imgaussfilt(tortuosity_mia_high_low.high_low.ch_one_power_transition_high_low,4),'linecolor','none');;
caxis([10 28]);
xlabel('High           Time from tortuosity switch (seconds)        Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Power (dB)");
title('MIA channel 1 (PFC)');
subplot(2,2,3)
contourf(power_times,power_freqs,imgaussfilt(tortuosity_control_high_low.high_low.ch_two_power_transition_high_low,4),'linecolor','none');
caxis([10 28]);
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Power (dB)");
title('Control channel 2 (HPC)');
subplot(2,2,4);
contourf(power_times,power_freqs,imgaussfilt(tortuosity_mia_high_low.high_low.ch_two_power_transition_high_low,4),'linecolor','none');
caxis([10 28]);
xlabel('High           Time from tortuosity switch (seconds)        Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Power (dB)");
title('MIA Channel 2 (HPC)');
sgtitle('Spectral power during high to low transitions');
colormap jet;

control_angle_mean_prechange_low_high = circ_mean(tortuosity_control_low_high.low_high.angle_prechange_session_low_high,[],2);
control_angle_mean_postchange_low_high = circ_mean(tortuosity_control_low_high.low_high.angle_postchange_session_low_high,[],2);

control_angle_mean_prechange_high_low = circ_mean(tortuosity_control_high_low.high_low.angle_prechange_session_high_low,[],2);
control_angle_mean_postchange_high_low = circ_mean(tortuosity_control_high_low.high_low.angle_postchange_session_high_low,[],2);



control_angle_mean_prechange_per_animal_low_high = circ_mean(tortuosity_control_low_high.low_high.angle_prechange_session_low_high,[],1);
control_angle_mean_postchange_per_animal_low_high = circ_mean(tortuosity_control_low_high.low_high.angle_postchange_session_low_high,[],1);

control_angle_mean_prechange_per_animal_high_low = circ_mean(tortuosity_control_high_low.high_low.angle_prechange_session_high_low,[],1);
control_angle_mean_postchange_per_animal_high_low = circ_mean(tortuosity_control_high_low.high_low.angle_postchange_session_high_low,[],1);

angle_per_animal_high = vertcat(control_angle_mean_postchange_per_animal_low_high',control_angle_mean_prechange_per_animal_high_low');
angle_per_animal_low = vertcat(control_angle_mean_prechange_per_animal_low_high',control_angle_mean_postchange_per_animal_high_low');




%% spectral granger
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/CTL_Tortuosity_spectral_granger_Results_v3.mat');
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/MIA_Tortuosity_spectral_granger_Results_v3.mat');

for z = 1:size(spectral_granger_results_ctl,2)

    hpc_to_pfc_low_high_ctl(1:80,1:1951,z) = spectral_granger_results_ctl(z).hpc_to_pfc_spec_low_to_high;
    pfc_to_hpc_low_high_ctl(1:80,1:1951,z) = spectral_granger_results_ctl(z).pfc_to_hpc_spec_low_to_high;
    hpc_to_pfc_high_low_ctl(1:80,1:1951,z) = spectral_granger_results_ctl(z).hpc_to_pfc_spec_high_to_low;
    pfc_to_hpc_high_low_ctl(1:80,1:1951,z) = spectral_granger_results_ctl(z).pfc_to_hpc_spec_high_to_low;
end

hpc_to_pfc_low_high_ctl = squeeze(nanmean(hpc_to_pfc_low_high_ctl,3));
pfc_to_hpc_low_high_ctl = squeeze(nanmean(pfc_to_hpc_low_high_ctl,3));
hpc_to_pfc_high_low_ctl = squeeze(nanmean(hpc_to_pfc_high_low_ctl,3));
pfc_to_hpc_high_low_ctl = squeeze(nanmean(pfc_to_hpc_high_low_ctl,3));

for z = 1:size(spectral_granger_results_mia,2)
    if isnan(sum(sum(spectral_granger_results_mia(z).hpc_to_pfc_spec_high_to_low)))
        remove_vector(z,1) = false;
    else
        remove_vector(z,1) = true;
    end
end
spectral_granger_results_mia = spectral_granger_results_mia(remove_vector);

for z = 1:size(spectral_granger_results_mia,2)

    hpc_to_pfc_low_high_mia(1:80,1:1951,z) = spectral_granger_results_mia(z).hpc_to_pfc_spec_low_to_high;
    pfc_to_hpc_low_high_mia(1:80,1:1951,z) = spectral_granger_results_mia(z).pfc_to_hpc_spec_low_to_high;
    hpc_to_pfc_high_low_mia(1:80,1:1951,z) = spectral_granger_results_mia(z).hpc_to_pfc_spec_high_to_low;
    pfc_to_hpc_high_low_mia(1:80,1:1951,z) = spectral_granger_results_mia(z).pfc_to_hpc_spec_high_to_low;
end

hpc_to_pfc_low_high_mia = squeeze(nanmean(hpc_to_pfc_low_high_mia,3));
pfc_to_hpc_low_high_mia = squeeze(nanmean(pfc_to_hpc_low_high_mia,3));
hpc_to_pfc_high_low_mia = squeeze(nanmean(hpc_to_pfc_high_low_mia,3));
pfc_to_hpc_high_low_mia = squeeze(nanmean(pfc_to_hpc_high_low_mia,3));


spec_frex = spectral_granger_results_ctl.frex;
spec_times = spectral_granger_results_ctl.time;
plasma_map = plasma(1500);
figure(4)
set(0, 'DefaultFigureRenderer', 'painters');
subplot(2,2,1)
imagesc(spec_times,spec_frex,flipud(hpc_to_pfc_low_high_ctl),[0 0.002]);
set(gca,'YDir','normal')
xlabel('Low           Time from tortuosity switch (seconds)         High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('Control, hippocampus -> pfc, low to high');
subplot(2,2,2);
imagesc(spec_times,spec_frex,flipud(pfc_to_hpc_low_high_ctl),[0 0.002]);
set(gca,'YDir','normal');
xlabel('Low           Time from tortuosity switch (seconds)        High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('Control, pfc -> hippocampus, low to high');
subplot(2,2,3)
imagesc(spec_times,spec_frex,flipud(hpc_to_pfc_high_low_ctl),[0 0.002]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('Control, hippocampus -> pfc, high to low');
subplot(2,2,4)
imagesc(spec_times,spec_frex,flipud(pfc_to_hpc_high_low_ctl),[0 0.002]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('Control, pfc -> hippocampus, high to low');
sgtitle('Spectral granger causality at tortuosity transitions, control animals');
colormap(plasma_map);



figure(5)
set(0, 'DefaultFigureRenderer', 'painters');
subplot(2,2,1)
imagesc(spec_times,spec_frex,flipud(hpc_to_pfc_low_high_mia),[0 0.002]);
set(gca,'YDir','normal')
xlabel('Low           Time from tortuosity switch (seconds)         High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('mia, hippocampus -> pfc, low to high');
subplot(2,2,2);
imagesc(spec_times,spec_frex,flipud(pfc_to_hpc_low_high_mia),[0 0.002]);
set(gca,'YDir','normal');
xlabel('Low           Time from tortuosity switch (seconds)        High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('mia, pfc -> hippocampus, low to high');
subplot(2,2,3)
imagesc(spec_times,spec_frex,flipud(hpc_to_pfc_high_low_mia),[0 0.002]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('mia, hippocampus -> pfc, high to low');
subplot(2,2,4)
imagesc(spec_times,spec_frex,flipud(pfc_to_hpc_high_low_mia),[0 0.002]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('mia, pfc -> hippocampus, high to low');
sgtitle('Spectral granger causality at tortuosity transitions, mia animals');
colormap(plasma_map);

for z = 1:size(spectral_granger_results_ctl,2)
    hpc_to_pfc_low_high_ctl_all(1:80,1:1951,z) = spectral_granger_results_ctl(z).hpc_to_pfc_spec_low_to_high;
    pfc_to_hpc_low_high_ctl_all(1:80,1:1951,z) = spectral_granger_results_ctl(z).pfc_to_hpc_spec_low_to_high;
    hpc_to_pfc_high_low_ctl_all(1:80,1:1951,z) = spectral_granger_results_ctl(z).hpc_to_pfc_spec_high_to_low;
    pfc_to_hpc_high_low_ctl_all(1:80,1:1951,z) = spectral_granger_results_ctl(z).pfc_to_hpc_spec_high_to_low;
end

ctl_h_to_p_lh_mean = nanmean(nanmean(hpc_to_pfc_low_high_ctl_all,3),1);
ctl_h_to_p_lh_sem = nanstd(nanstd(hpc_to_pfc_low_high_ctl_all,0,3),0,1)/sqrt(size(spectral_granger_results_ctl,2));

ctl_p_to_h_lh_mean = nanmean(nanmean(pfc_to_hpc_low_high_ctl_all,3),1);
ctl_p_to_h_lh_sem = nanstd(nanstd(pfc_to_hpc_low_high_ctl_all,0,3),0,1)/sqrt(size(spectral_granger_results_ctl,2));

ctl_h_to_p_hl_mean = nanmean(nanmean(hpc_to_pfc_high_low_ctl_all,3),1);
ctl_h_to_p_hl_sem = nanstd(nanstd(hpc_to_pfc_high_low_ctl_all,0,3),0,1)/sqrt(size(spectral_granger_results_ctl,2));

ctl_p_to_h_hl_mean = nanmean(nanmean(pfc_to_hpc_high_low_ctl_all,3),1);
ctl_p_to_h_hl_sem = nanstd(nanstd(pfc_to_hpc_high_low_ctl_all,0,3),0,1)/sqrt(size(spectral_granger_results_ctl,2));

for z = 1:size(spectral_granger_results_mia,2)
    hpc_to_pfc_low_high_mia_all(1:80,1:1951,z) = spectral_granger_results_mia(z).hpc_to_pfc_spec_low_to_high;
    pfc_to_hpc_low_high_mia_all(1:80,1:1951,z) = spectral_granger_results_mia(z).pfc_to_hpc_spec_low_to_high;
    hpc_to_pfc_high_low_mia_all(1:80,1:1951,z) = spectral_granger_results_mia(z).hpc_to_pfc_spec_high_to_low;
    pfc_to_hpc_high_low_mia_all(1:80,1:1951,z) = spectral_granger_results_mia(z).pfc_to_hpc_spec_high_to_low;
end


mia_h_to_p_lh_mean = nanmean(nanmean(hpc_to_pfc_low_high_mia,3),1);
mia_h_to_p_lh_sem = nanstd(nanstd(hpc_to_pfc_low_high_mia,0,3),0,1)/sqrt(size(spectral_granger_results_mia,2));

mia_p_to_h_lh_mean = nanmean(nanmean(pfc_to_hpc_low_high_mia,3),1);
mia_p_to_h_lh_sem = nanstd(nanstd(pfc_to_hpc_low_high_mia,0,3),0,1)/sqrt(size(spectral_granger_results_mia,2));

mia_h_to_p_hl_mean = nanmean(nanmean(hpc_to_pfc_high_low_mia,3),1);
mia_h_to_p_hl_sem = nanstd(nanstd(hpc_to_pfc_high_low_mia,0,3),0,1)/sqrt(size(spectral_granger_results_mia,2));

mia_p_to_h_hl_mean = nanmean(nanmean(pfc_to_hpc_high_low_mia,3),1);
mia_p_to_h_hl_sem = nanstd(nanstd(pfc_to_hpc_high_low_mia,0,3),0,1)/sqrt(size(spectral_granger_results_mia,2));

figure(6)
shadedErrorBar(spec_times(1:1951),ctl_h_to_p_lh_mean,ctl_h_to_p_lh_sem);
hold on
shadedErrorBar(spec_times(1:1951),mia_h_to_p_lh_mean,mia_h_to_p_lh_sem,'lineProps','r-');
xline(0,'k--','LineWidth',2);
ylabel('Granger causality');
xlabel('Low         Time from transition (s)        High');
title('hpc -> pfc low to high');

figure(7)
shadedErrorBar(spec_times(1:1951),ctl_p_to_h_lh_mean,ctl_p_to_h_lh_sem);
hold on
shadedErrorBar(spec_times(1:1951),mia_p_to_h_lh_mean,mia_p_to_h_lh_sem,'lineProps','r-');
xline(0,'k--','LineWidth',2);
ylabel('Granger causality');
xlabel('Low         Time from transition (s)        High');
title('pfc -> hpc low to high');

figure(8)
shadedErrorBar(spec_times(1:1951),ctl_h_to_p_hl_mean,ctl_h_to_p_hl_sem);
hold on
shadedErrorBar(spec_times(1:1951),mia_h_to_p_hl_mean,mia_h_to_p_hl_sem,'lineProps','r-');
xline(0,'k--','LineWidth',2);
ylabel('Granger causality');
xlabel('High         Time from transition (s)        Low');
title('hpc -> pfc high to low');

figure(9)
shadedErrorBar(spec_times(1:1951),ctl_p_to_h_hl_mean,ctl_p_to_h_hl_sem);
hold on
shadedErrorBar(spec_times(1:1951),mia_p_to_h_hl_mean,mia_p_to_h_hl_sem,'lineProps','r-');
xline(0,'k--','LineWidth',2);
ylabel('Granger causality');
xlabel('High        Time from transition (s)        Low');
title('pfc -> hpc high to low');




ctl_h_to_p_lh = nanmean(hpc_to_pfc_low_high_ctl_all,3);
ctl_p_to_h_lh = nanmean(pfc_to_hpc_low_high_ctl_all,3);
ctl_h_to_p_hl = nanmean(hpc_to_pfc_high_low_ctl_all,3);
ctl_p_to_h_hl = nanmean(pfc_to_hpc_high_low_ctl_all,3);


figure(10)
set(0, 'DefaultFigureRenderer', 'painters');
subplot(2,2,1)
imagesc(spec_times,spec_frex,flipud(ctl_h_to_p_lh),[0 0.004]);
set(gca,'YDir','normal')
xlabel('Low           Time from tortuosity switch (seconds)         High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('Control, hippocampus -> pfc, low to high');
subplot(2,2,2);
imagesc(spec_times,spec_frex,flipud(ctl_p_to_h_lh),[0 0.004]);
set(gca,'YDir','normal');
xlabel('Low           Time from tortuosity switch (seconds)        High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('Control, pfc -> hippocampus, low to high');
subplot(2,2,3)
imagesc(spec_times,spec_frex,flipud(ctl_h_to_p_hl),[0 0.004]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('Control, hippocampus -> pfc, high to low');
subplot(2,2,4)
imagesc(spec_times,spec_frex,flipud(ctl_p_to_h_hl),[0 0.004]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('Control, pfc -> hippocampus, high to low');
sgtitle('Spectral granger causality at tortuosity transitions, control animals');
colormap(plasma_map);

mia_h_to_p_lh = nanmean(hpc_to_pfc_low_high_mia,3);
mia_p_to_h_lh = nanmean(pfc_to_hpc_low_high_mia,3);
mia_h_to_p_hl = nanmean(hpc_to_pfc_high_low_mia,3);
mia_p_to_h_hl = nanmean(pfc_to_hpc_high_low_mia,3);


figure(11)
set(0, 'DefaultFigureRenderer', 'painters');
subplot(2,2,1)
imagesc(spec_times,spec_frex,flipud(mia_h_to_p_lh),[0 0.004]);
set(gca,'YDir','normal')
xlabel('Low           Time from tortuosity switch (seconds)         High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('MIA, hippocampus -> pfc, low to high');
subplot(2,2,2);
imagesc(spec_times,spec_frex,flipud(mia_p_to_h_lh),[0 0.004]);
set(gca,'YDir','normal');
xlabel('Low           Time from tortuosity switch (seconds)        High');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('MIA, pfc -> hippocampus, low to high');
subplot(2,2,3)
imagesc(spec_times,spec_frex,flipud(mia_h_to_p_hl),[0 0.004]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('MIA, hippocampus -> pfc, high to low');
subplot(2,2,4)
imagesc(spec_times,spec_frex,flipud(mia_p_to_h_hl),[0 0.004]);
set(gca,'YDir','normal')
xlabel('High           Time from tortuosity switch (seconds)         Low');
ylabel('Frequency (Hz)');
xline(0,'k--','LineWidth',4);
hcb = colorbar;
ylabel(hcb,"Granger causality");
title('MIA, pfc -> hippocampus, high to low');
sgtitle('Spectral granger causality at tortuosity transitions, MIA animals');
colormap(plasma_map);


%% By animal

%%%% coherence
%%%% some of these coherence values are whack low, so we need to weed them
%%%% out

for i = 1:size(tortuosity_control_high_low.high_low.coherence_transition_high_low_by_animal,3)
    mean_mat = nanmean(tortuosity_control_high_low.high_low.coherence_transition_high_low_by_animal,'all');

    hi_low = tortuosity_control_high_low.high_low.coherence_transition_high_low_by_animal(:,:,i);
    hi_low(hi_low > 2*mean_mat) = mean_mat;
    hi_low(hi_low < -2*mean_mat) = mean_mat;
    tortuosity_control_high_low.high_low.coherence_transition_high_low_by_animal(:,:,i) = hi_low;

    low_hi = tortuosity_control_low_high.low_high.coherence_transition_low_high_by_animal(:,:,i);
    low_hi(low_hi > 2*mean_mat) = mean_mat;
    low_hi(low_hi < -2*mean_mat) = mean_mat;
    tortuosity_control_low_high.low_high.coherence_transition_low_high_by_animal(:,:,i) = low_hi;
end

high_low_control_coherence_prechange = squeeze(nanmean(tortuosity_control_high_low.high_low.coherence_transition_high_low_by_animal,1));
low_high_control_coherence_prechange = squeeze(nanmean(tortuosity_control_low_high.low_high.coherence_transition_low_high_by_animal,1));


for z = 1:size(tortuosity_control_high_low.high_low.coherence_transition_high_low_by_animal,3)
    hi_lo_pre(z,1) = nanmean(high_low_coherence_prechange(500:1000,z));
    hi_lo_post(z,1) = nanmean(high_low_coherence_prechange(1001:1501,z));
end

for z = 1:size(tortuosity_control_low_high.high_low.coherence_transition_high_low_by_animal,3)
    hi_lo_pre(z,1) = nanmean(high_low_coherence_prechange(500:1000,z));
    hi_lo_post(z,1) = nanmean(high_low_coherence_prechange(1001:1501,z));
end

%%%%%%%%%%%%%%%%%%%%
%% COHEN ANALYSES %%
%%%%%%%%%%%%%%%%%%%%

load('/Volumes/1TB ALPHA/MIA_Coherence_Results/CTL_master_cohen_results_v3.mat')
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/MIA_master_cohen_results_v3.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ERP ANALYSIS %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONTROL
%%%% LOW TO HIGH TRANSITIONS
for i = 1:size(erp_results_high_low_CTL,2)
    hl_erp_ch_one(i,:) = erp_results_high_low_CTL(i).ERP_ch_one;
    hl_erp_ch_two(i,:) = erp_results_high_low_CTL(i).erp_ch_two;

    lh_erp_ch_one(i,:) = erp_results_low_high_CTL(i).ERP_ch_one;
    lh_erp_ch_two(i,:) = erp_results_low_high_CTL(i).erp_ch_two;
end
times = linspace(-4,4,2001);

for i = 1:size(erp_results_low_high_CTL,2)
    total_low_high_channel_one(i,:,:) = erp_results_low_high_CTL(i).total_power_ch_one;
    total_low_high_channel_two(i,:,:) = erp_results_low_high_CTL(i).total_power_ch_two;

    non_phase_locked_low_high_channel_one(i,:,:) = erp_results_low_high_CTL(i).non_phase_locked_power_ch_one;
    non_phase_locked_low_high_channel_two(i,:,:) = erp_results_low_high_CTL(i).non_phase_locked_power_ch_two;

    erp_power_low_high_channel_one(i,:,:) = erp_results_low_high_CTL(i).ERP_power_ch_one;
    erp_power_low_high_channel_two(i,:,:) = erp_results_low_high_CTL(i).ERP_power_ch_two;

    itpc_low_high_channel_one(i,:,:) = erp_results_low_high_CTL(i).ITPC_ch_one;
    itpc_low_high_channel_two(i,:,:) = erp_results_low_high_CTL(i).ITPC_ch_two;
end

%%%% GET CONTROL DATA

for i = 1:size(erp_results_low_high_CTL,2)
    total_channel_one_control(i,:,:) = erp_results_low_high_CTL(i).total_power_ch_one_control;
    total_channel_two_control(i,:,:) = erp_results_low_high_CTL(i).total_power_ch_two_control;

    non_phase_locked_channel_one_control(i,:,:) = erp_results_low_high_CTL(i).non_phase_locked_power_ch_one_control;
    non_phase_locked_channel_two_control(i,:,:) = erp_results_low_high_CTL(i).non_phase_locked_power_ch_two_control;

    erp_power_channel_one_control(i,:,:) = erp_results_low_high_CTL(i).ERP_power_ch_one_control;
    erp_power_channel_two_control(i,:,:) = erp_results_low_high_CTL(i).ERP_power_ch_two_control;

    itpc_channel_one_control(i,:,:) = erp_results_low_high_CTL(i).ITPC_ch_one_control;
    itpc_channel_two_control(i,:,:) = erp_results_low_high_CTL(i).ITPC_ch_two_control;
end

analysis_labels = {'Total';'Non-phase-locked';'ERP power';'ITPC'};
frequencies = erp_results_low_high_CTL(1).freqs;

erp_data_channel_one_low_high_transitions(1,:,:) = squeeze(nanmean(total_low_high_channel_one));
erp_data_channel_one_low_high_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_low_high_channel_one));
erp_data_channel_one_low_high_transitions(3,:,:) = squeeze(nanmean(erp_power_low_high_channel_one));
erp_data_channel_one_low_high_transitions(4,:,:) = squeeze(nanmean(itpc_low_high_channel_one));


erp_control(1,:,:) = squeeze(nanmean(total_channel_one_control));
erp_control(2,:,:) = squeeze(nanmean(non_phase_locked_channel_one_control));
erp_control(3,:,:) = squeeze(nanmean(erp_power_channel_one_control));
erp_control(4,:,:) = squeeze(nanmean(itpc_channel_one_control));

% color limits
clims = [ -1 1; -1 1; -3 3; 0.1 .2 ];
erp = nanmean(lh_erp_ch_one);
% scale ERP for plotting
erpt = (erp-min(erp))./max(erp-min(erp));
erpt = erpt*(frequencies(end)-frequencies(1))+frequencies(1);

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_data_channel_one_low_high_transitions(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('Low         Time from switch (ms)         High')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel one, low to high transitions')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_control(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('Low         Time from switch (ms)         High')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel one, control segments')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

%%% CHANNEL TWO
erp_data_channel_two_low_high_transitions(1,:,:) = squeeze(nanmean(total_low_high_channel_two));
erp_data_channel_two_low_high_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_low_high_channel_two));
erp_data_channel_two_low_high_transitions(3,:,:) = squeeze(nanmean(erp_power_low_high_channel_two));
erp_data_channel_two_low_high_transitions(4,:,:) = squeeze(nanmean(itpc_low_high_channel_two));

erp_control_two(1,:,:) = squeeze(nanmean(total_channel_two_control));
erp_control_two(2,:,:) = squeeze(nanmean(non_phase_locked_channel_two_control));
erp_control_two(3,:,:) = squeeze(nanmean(erp_power_channel_two_control));
erp_control_two(4,:,:) = squeeze(nanmean(itpc_channel_two_control));

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_data_channel_two_low_high_transitions(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('Low         Time from switch (ms)         High')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel two, low to high transitions')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');


figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_control_two(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('Low         Time from switch (ms)         High')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel two, control segments ')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

%%%%% HIGH TO LOW TRANSITIONS
for i = 1:size(erp_results_high_low_CTL,2)
    total_high_low_channel_one(i,:,:) = erp_results_high_low_CTL(i).total_power_ch_one;
    total_high_low_channel_two(i,:,:) = erp_results_high_low_CTL(i).total_power_ch_two;

    non_phase_locked_high_low_channel_one(i,:,:) = erp_results_high_low_CTL(i).non_phase_locked_power_ch_one;
    non_phase_locked_high_low_channel_two(i,:,:) = erp_results_high_low_CTL(i).non_phase_locked_power_ch_two;

    erp_power_high_low_channel_one(i,:,:) = erp_results_high_low_CTL(i).ERP_power_ch_one;
    erp_power_high_low_channel_two(i,:,:) = erp_results_high_low_CTL(i).ERP_power_ch_two;

    itpc_high_low_channel_one(i,:,:) = erp_results_high_low_CTL(i).ITPC_ch_one;
    itpc_high_low_channel_two(i,:,:) = erp_results_high_low_CTL(i).ITPC_ch_two;
end

for i = 1:size(erp_results_high_low_CTL,2)
    total_channel_one_control(i,:,:) = erp_results_high_low_CTL(i).total_power_ch_one_control;
    total_channel_two_control(i,:,:) = erp_results_high_low_CTL(i).total_power_ch_two_control;

    non_phase_locked_channel_one_control(i,:,:) = erp_results_high_low_CTL(i).non_phase_locked_power_ch_one_control;
    non_phase_locked_channel_two_control(i,:,:) = erp_results_high_low_CTL(i).non_phase_locked_power_ch_two_control;

    erp_power_channel_one_control(i,:,:) = erp_results_high_low_CTL(i).ERP_power_ch_one_control;
    erp_power_channel_two_control(i,:,:) = erp_results_high_low_CTL(i).ERP_power_ch_two_control;

    itpc_channel_one_control(i,:,:) = erp_results_high_low_CTL(i).ITPC_ch_one_control;
    itpc_channel_two_control(i,:,:) = erp_results_high_low_CTL(i).ITPC_ch_two_control;
end

erp_control_one(1,:,:) = squeeze(nanmean(total_channel_one_control));
erp_control_one(2,:,:) = squeeze(nanmean(non_phase_locked_channel_one_control));
erp_control_one(3,:,:) = squeeze(nanmean(erp_power_channel_one_control));
erp_control_one(4,:,:) = squeeze(nanmean(itpc_channel_one_control));

erp_control_two(1,:,:) = squeeze(nanmean(total_channel_two_control));
erp_control_two(2,:,:) = squeeze(nanmean(non_phase_locked_channel_two_control));
erp_control_two(3,:,:) = squeeze(nanmean(erp_power_channel_two_control));
erp_control_two(4,:,:) = squeeze(nanmean(itpc_channel_two_control));

analysis_labels = {'Total';'Non-phase-locked';'ERP power';'ITPC'};
frequencies = erp_results_high_low_CTL(1).freqs;

erp_data_channel_one_high_low_transitions(1,:,:) = squeeze(nanmean(total_high_low_channel_one));
erp_data_channel_one_high_low_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_high_low_channel_one));
erp_data_channel_one_high_low_transitions(3,:,:) = squeeze(nanmean(erp_power_high_low_channel_one));
erp_data_channel_one_high_low_transitions(4,:,:) = squeeze(nanmean(itpc_high_low_channel_one));

erp_data_channel_two_high_low_transitions(1,:,:) = squeeze(nanmean(total_high_low_channel_two));
erp_data_channel_two_high_low_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_high_low_channel_two));
erp_data_channel_two_high_low_transitions(3,:,:) = squeeze(nanmean(erp_power_high_low_channel_two));
erp_data_channel_two_high_low_transitions(4,:,:) = squeeze(nanmean(itpc_high_low_channel_two));

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_data_channel_one_high_low_transitions(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('High         Time from switch (ms)         Low')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel one, high to low transitions')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_data_channel_two_high_low_transitions(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4,'yscale','log')
    xline(0,'k--');
    xlabel('High        Time from switch (ms)         Low')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel two, high to low transitions')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_control_one(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('High         Time from switch (ms)         Low')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel one, control segments')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_control_two(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4,'yscale','log')
    xline(0,'k--');
    xlabel('High        Time from switch (ms)         Low')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel two, control segments')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

erp_theta_range = 6:11;
for z = 1:size(total_low_high_channel_one,1)
    total_low_high_ch_one_theta_prechange(z,1) = mean(mean(total_low_high_channel_one(z,erp_theta_range,1:1000)));
    total_low_high_ch_one_theta_postchange(z,1) = mean(mean(total_low_high_channel_one(z,erp_theta_range,1001:2001)));
end

%% MIA
clear  total_low_high_channel_one  total_low_high_channel_two non_phase_locked_low_high_channel_one non_phase_locked_low_high_channel_two...
    erp_power_low_high_channel_one erp_power_low_high_channel_two itpc_low_high_channel_one itpc_low_high_channel_two
for i = 1:size(erp_results_high_low_mia,2)
    hl_erp_ch_one(i,:) = erp_results_high_low_mia(i).ERP_ch_one;
    hl_erp_ch_two(i,:) = erp_results_high_low_mia(i).erp_ch_two;

    lh_erp_ch_one(i,:) = erp_results_low_high_mia(i).ERP_ch_one;
    lh_erp_ch_two(i,:) = erp_results_low_high_mia(i).erp_ch_two;
end
times = linspace(-4,4,2001);

for i = 1:size(erp_results_low_high_mia,2)
    total_low_high_channel_one(i,:,:) = erp_results_low_high_mia(i).total_power_ch_one;
    total_low_high_channel_two(i,:,:) = erp_results_low_high_mia(i).total_power_ch_two;

    non_phase_locked_low_high_channel_one(i,:,:) = erp_results_low_high_mia(i).non_phase_locked_power_ch_one;
    non_phase_locked_low_high_channel_two(i,:,:) = erp_results_low_high_mia(i).non_phase_locked_power_ch_two;

    erp_power_low_high_channel_one(i,:,:) = erp_results_low_high_mia(i).ERP_power_ch_one;
    erp_power_low_high_channel_two(i,:,:) = erp_results_low_high_mia(i).ERP_power_ch_two;

    itpc_low_high_channel_one(i,:,:) = erp_results_low_high_mia(i).ITPC_ch_one;
    itpc_low_high_channel_two(i,:,:) = erp_results_low_high_mia(i).ITPC_ch_two;
end

%%%% GET MIA DATA
clear total_channel_one_control total_channel_two_control  non_phase_locked_channel_one_control  non_phase_locked_channel_two_control...
     erp_power_channel_one_control  erp_power_channel_two_control itpc_channel_one_control itpc_channel_two_control...
     erp_data_channel_one_low_high_transitions erp_data_channel_two_low_high_transitions...
     erp_data_channel_two_high_low_transitions erp_data_channel_two_high_low_transitions...
     erp_control erp_control_two erp_data_channel_two_low_high_transitions

for i = 1:size(erp_results_low_high_mia,2)
    total_channel_one_control(i,:,:) = erp_results_low_high_mia(i).total_power_ch_one_control;
    total_channel_two_control(i,:,:) = erp_results_low_high_mia(i).total_power_ch_two_control;

    non_phase_locked_channel_one_control(i,:,:) = erp_results_low_high_mia(i).non_phase_locked_power_ch_one_control;
    non_phase_locked_channel_two_control(i,:,:) = erp_results_low_high_mia(i).non_phase_locked_power_ch_two_control;

    erp_power_channel_one_control(i,:,:) = erp_results_low_high_mia(i).ERP_power_ch_one_control;
    erp_power_channel_two_control(i,:,:) = erp_results_low_high_mia(i).ERP_power_ch_two_control;

    itpc_channel_one_control(i,:,:) = erp_results_low_high_mia(i).ITPC_ch_one_control;
    itpc_channel_two_control(i,:,:) = erp_results_low_high_mia(i).ITPC_ch_two_control;
end

analysis_labels = {'Total';'Non-phase-locked';'ERP power';'ITPC'};
frequencies = erp_results_low_high_mia(1).freqs;

erp_data_channel_one_low_high_transitions(1,:,:) = squeeze(nanmean(total_low_high_channel_one));
erp_data_channel_one_low_high_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_low_high_channel_one));
erp_data_channel_one_low_high_transitions(3,:,:) = squeeze(nanmean(erp_power_low_high_channel_one));
erp_data_channel_one_low_high_transitions(4,:,:) = squeeze(nanmean(itpc_low_high_channel_one));


erp_control(1,:,:) = squeeze(nanmean(total_channel_one_control));
erp_control(2,:,:) = squeeze(nanmean(non_phase_locked_channel_one_control));
erp_control(3,:,:) = squeeze(nanmean(erp_power_channel_one_control));
erp_control(4,:,:) = squeeze(nanmean(itpc_channel_one_control));

% color limits
clims = [ -1 1; -1 1; -3 3; 0.1 .2 ];
erp = nanmean(lh_erp_ch_one);
% scale ERP for plotting
erpt = (erp-min(erp))./max(erp-min(erp));
erpt = erpt*(frequencies(end)-frequencies(1))+frequencies(1);

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_data_channel_one_low_high_transitions(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('Low         Time from switch (ms)         High')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel one, low to high transitions')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_control(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('Low         Time from switch (ms)         High')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel one, control segments')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

%%% CHANNEL TWO
erp_data_channel_two_low_high_transitions(1,:,:) = squeeze(nanmean(total_low_high_channel_two));
erp_data_channel_two_low_high_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_low_high_channel_two));
erp_data_channel_two_low_high_transitions(3,:,:) = squeeze(nanmean(erp_power_low_high_channel_two));
erp_data_channel_two_low_high_transitions(4,:,:) = squeeze(nanmean(itpc_low_high_channel_two));

erp_control_two(1,:,:) = squeeze(nanmean(total_channel_two_control));
erp_control_two(2,:,:) = squeeze(nanmean(non_phase_locked_channel_two_control));
erp_control_two(3,:,:) = squeeze(nanmean(erp_power_channel_two_control));
erp_control_two(4,:,:) = squeeze(nanmean(itpc_channel_two_control));

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_data_channel_two_low_high_transitions(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('Low         Time from switch (ms)         High')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel two, low to high transitions')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');


figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_control_two(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('Low         Time from switch (ms)         High')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel two, control segments ')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

%%%%% HIGH TO LOW TRANSITIONS
for i = 1:size(erp_results_high_low_mia,2)
    total_high_low_channel_one(i,:,:) = erp_results_high_low_mia(i).total_power_ch_one;
    total_high_low_channel_two(i,:,:) = erp_results_high_low_mia(i).total_power_ch_two;

    non_phase_locked_high_low_channel_one(i,:,:) = erp_results_high_low_mia(i).non_phase_locked_power_ch_one;
    non_phase_locked_high_low_channel_two(i,:,:) = erp_results_high_low_mia(i).non_phase_locked_power_ch_two;

    erp_power_high_low_channel_one(i,:,:) = erp_results_high_low_mia(i).ERP_power_ch_one;
    erp_power_high_low_channel_two(i,:,:) = erp_results_high_low_mia(i).ERP_power_ch_two;

    itpc_high_low_channel_one(i,:,:) = erp_results_high_low_mia(i).ITPC_ch_one;
    itpc_high_low_channel_two(i,:,:) = erp_results_high_low_mia(i).ITPC_ch_two;
end

for i = 1:size(erp_results_high_low_mia,2)
    total_channel_one_control(i,:,:) = erp_results_high_low_mia(i).total_power_ch_one_control;
    total_channel_two_control(i,:,:) = erp_results_high_low_mia(i).total_power_ch_two_control;

    non_phase_locked_channel_one_control(i,:,:) = erp_results_high_low_mia(i).non_phase_locked_power_ch_one_control;
    non_phase_locked_channel_two_control(i,:,:) = erp_results_high_low_mia(i).non_phase_locked_power_ch_two_control;

    erp_power_channel_one_control(i,:,:) = erp_results_high_low_mia(i).ERP_power_ch_one_control;
    erp_power_channel_two_control(i,:,:) = erp_results_high_low_mia(i).ERP_power_ch_two_control;

    itpc_channel_one_control(i,:,:) = erp_results_high_low_mia(i).ITPC_ch_one_control;
    itpc_channel_two_control(i,:,:) = erp_results_high_low_mia(i).ITPC_ch_two_control;
end

erp_control_one(1,:,:) = squeeze(nanmean(total_channel_one_control));
erp_control_one(2,:,:) = squeeze(nanmean(non_phase_locked_channel_one_control));
erp_control_one(3,:,:) = squeeze(nanmean(erp_power_channel_one_control));
erp_control_one(4,:,:) = squeeze(nanmean(itpc_channel_one_control));

erp_control_two(1,:,:) = squeeze(nanmean(total_channel_two_control));
erp_control_two(2,:,:) = squeeze(nanmean(non_phase_locked_channel_two_control));
erp_control_two(3,:,:) = squeeze(nanmean(erp_power_channel_two_control));
erp_control_two(4,:,:) = squeeze(nanmean(itpc_channel_two_control));

analysis_labels = {'Total';'Non-phase-locked';'ERP power';'ITPC'};
frequencies = erp_results_high_low_mia(1).freqs;

erp_data_channel_one_high_low_transitions(1,:,:) = squeeze(nanmean(total_high_low_channel_one));
erp_data_channel_one_high_low_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_high_low_channel_one));
erp_data_channel_one_high_low_transitions(3,:,:) = squeeze(nanmean(erp_power_high_low_channel_one));
erp_data_channel_one_high_low_transitions(4,:,:) = squeeze(nanmean(itpc_high_low_channel_one));

erp_data_channel_two_high_low_transitions(1,:,:) = squeeze(nanmean(total_high_low_channel_two));
erp_data_channel_two_high_low_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_high_low_channel_two));
erp_data_channel_two_high_low_transitions(3,:,:) = squeeze(nanmean(erp_power_high_low_channel_two));
erp_data_channel_two_high_low_transitions(4,:,:) = squeeze(nanmean(itpc_high_low_channel_two));

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_data_channel_one_high_low_transitions(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('High         Time from switch (ms)         Low')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel one, high to low transitions')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_data_channel_two_high_low_transitions(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4,'yscale','log')
    xline(0,'k--');
    xlabel('High        Time from switch (ms)         Low')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel two, high to low transitions')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_control_one(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
    xline(0,'k--');
    xlabel('High         Time from switch (ms)         Low')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel one, control segments')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

figure("Position",[10 10 2000 600],"Visible","on")
for i=1:4
    subplot(2,3,i)
    contourf(times,frequencies,squeeze(erp_control_two(i,:,:)),40,'linecolor','none')
    set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4,'yscale','log')
    xline(0,'k--');
    xlabel('High        Time from switch (ms)         Low')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    sgtitle('Channel two, control segments')
    if i == 3
        c = colorbar();
ylabel(c,'Baseline-subtracted power (dB)');
    end
   
    colormap("jet")
end
        c = colorbar();
ylabel(c,'Baseline-subtracted Phase clustering');

erp_theta_range = 6:11;
for z = 1:size(total_low_high_channel_one,1)
    total_low_high_ch_one_theta_prechange(z,1) = mean(mean(total_low_high_channel_one(z,erp_theta_range,1:1000)));
    total_low_high_ch_one_theta_postchange(z,1) = mean(mean(total_low_high_channel_one(z,erp_theta_range,1001:2001)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% COHERENCE ANALYSES %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:size(cohen_cohere_results_low_high_CTL,2)
    bs_low_high_coherence(ii,:,:) = cohen_cohere_results_low_high_CTL(ii).baseline_subtracted_imaginary_coherence;
    bs_high_low_coherence(ii,:,:) = cohen_cohere_results_high_low_CTL(ii).baseline_subtracted_imaginary_coherence;
    low_high_coherence(ii,:,:) = cohen_cohere_results_low_high_CTL(ii).coherence;
    high_low_coherence(ii,:,:) = cohen_cohere_results_high_low_CTL(ii).coherence;
    low_high_imag_coherence(ii,:,:) = cohen_cohere_results_low_high_CTL(ii).imaginary_coherence;
    high_low_imag_coherence(ii,:,:) = cohen_cohere_results_high_low_CTL(ii).imaginary_coherence;
    control_coherence_imag_lh(ii,:,:) = cohen_cohere_results_low_high_CTL(ii).imaginary_coherence_control;
    control_coherence_imag_hl(ii,:,:) = cohen_cohere_results_high_low_CTL(ii).imaginary_coherence_control;
    control_coherence_lh(ii,:,:) = cohen_cohere_results_low_high_CTL(ii).coherence_control;
    control_coherence_hl(ii,:,:) = cohen_cohere_results_high_low_CTL(ii).coherence_control;
end

bsub_low_high_coherence = squeeze(mean(bs_low_high_coherence,1));
bsub_high_low_coherence = squeeze(mean(bs_high_low_coherence,1));
lh_cohere = squeeze(mean(low_high_coherence,1));
hl_cohere = squeeze(mean(high_low_coherence,1));
lh_cohere_control = squeeze(mean(control_coherence_lh,1));
hl_cohere_control = squeeze(mean(control_coherence_hl,1));
lh_imag_cohere = squeeze(mean(low_high_imag_coherence,1));
hl_imag_cohere = squeeze(mean(high_low_imag_coherence,1));
lh_imag_cohere_control = squeeze(mean(control_coherence_imag_lh,1));
hl_imag_cohere_control = squeeze(mean(control_coherence_imag_hl,1));

cohere_times = linspace(-4,4,161);
cohere_freqs = cohen_cohere_results_low_high_CTL(1).freqs;
rdylbu = brewermap(1024,'-RdYlBu');
ylorrd = brewermap(1024,'YlOrRd');

%%% baseline subtracted coherence
figure("Position",[10 10 2000 600],"Visible","on")
subplot(1,2,1)
sgtitle('Baseline-subtracted coherence at transitions')
contourf(cohere_times,cohere_freqs,bsub_low_high_coherence,40,'linecolor','none')
set(gca,'clim',[0 0.2],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xline(0,'k--');
xlabel('Low           Time from regime change (s)           High');
ylabel('Frequency (Hz)')
title('Low to High transitions')
colormap("jet")
subplot(1,2,2)
contourf(cohere_times,cohere_freqs,bsub_high_low_coherence,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('High           Time from regime change (s)           Low');
ylabel('Frequency (Hz)')
title('High to Low transitions')
c = colorbar();
ylabel(c,'Baseline subtracted coherence')
colormap("jet")

%%% raw imaginary coherence
figure("Position",[10 10 2000 600],"Visible","on")
sgtitle('Coherence at transitions')
subplot(1,2,1)
contourf(cohere_times,cohere_freqs,lh_imag_cohere,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('Low           Time from regime change (s)           High');
ylabel('Frequency (Hz)')
colormap("jet")
subplot(1,2,2)
contourf(cohere_times,cohere_freqs,hl_imag_cohere,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('High           Time from regime change (s)           Low');
ylabel('Frequency (Hz)')
c = colorbar();
ylabel(c,'Coherence');
colormap("jet")

%%% raw imaginary coherence
figure("Position",[10 10 2000 600],"Visible","on")
sgtitle('Coherence at transitions')
subplot(1,2,1)
contourf(cohere_times,cohere_freqs,lh_imag_cohere_control,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('Low           Time from regime change (s)           High');
ylabel('Frequency (Hz)')
colormap("jet")
subplot(1,2,2)
contourf(cohere_times,cohere_freqs,hl_imag_cohere_control,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('High           Time from regime change (s)           Low');
ylabel('Frequency (Hz)')
c = colorbar();
ylabel(c,'Coherence');
colormap("jet")



time_around_transition = 71:91;
one_to_two_sec_prior = 1:21;
theta_coherence_freqs = 27:52;

for z = 1:size(bs_low_high_coherence,1)
    ctl_theta_coherence_low_high_around_transition(z,:) = mean(squeeze(mean(bs_low_high_coherence(z,theta_coherence_freqs,time_around_transition))));
    ctl_theta_coherence_low_high_four_to_three(z,:) = mean(squeeze(mean(bs_low_high_coherence(z,theta_coherence_freqs,one_to_two_sec_prior))));

    ctl_theta_coherence_high_low_around_transition(z,:) = mean(squeeze(mean(bs_high_low_coherence(z,theta_coherence_freqs,time_around_transition))));
    ctl_theta_coherence_high_low_four_to_three(z,:) = mean(squeeze(mean(bs_high_low_coherence(z,theta_coherence_freqs,one_to_two_sec_prior))));
end

low_high_control_theta_coherence = vertcat(ctl_theta_coherence_low_high_four_to_three,ctl_theta_coherence_low_high_around_transition);
high_low_control_theta_coherence = vertcat(ctl_theta_coherence_high_low_four_to_three,ctl_theta_coherence_high_low_around_transition);
type(1:42,1) = {'Pre-transition'}; type(43:42*2,1) = {'Peri-Transition'};
axis(1:42,1) = 1:42; axis(43:42*2) = 1:42;

figure()
g = gramm('x',type,'y',low_high_control_theta_coherence,'group',axis);
g.geom_line();
g.set_order_options("x",0)
g.update('x',type,'y',low_high_control_theta_coherence,'color',type)
g.geom_point();
g.set_names('x','','y','Baseline-subtracted coherence');
g.set_title('Low to high transitions');
g.axe_property('YLim',[0 0.35]);
g.draw();

figure()
g = gramm('x',type,'y',high_low_control_theta_coherence,'group',axis);
g.geom_line();
g.set_order_options("x",0)
g.update('x',type,'y',high_low_control_theta_coherence,'color',type)
g.geom_point();
g.set_names('x','','y','Baseline-subtracted coherence');
g.set_title('Low to high transitions');
g.axe_property('YLim',[0 0.35]);
g.draw();


[ctl_lh_coherence_p,~,ctl_lh_coherence_stat] = signrank(ctl_theta_coherence_low_high_four_to_three,ctl_theta_coherence_low_high_around_transition);

[ctl_hl_coherence_p,~,ctl_hl_coherence_stat] = signrank(ctl_theta_coherence_high_low_four_to_three,ctl_theta_coherence_high_low_around_transition);




%%%% MIA

for ii = 1:size(cohen_cohere_results_low_high_mia,2)
    bs_low_high_coherence(ii,:,:) = cohen_cohere_results_low_high_mia(ii).baseline_subtracted_imaginary_coherence;
    bs_high_low_coherence(ii,:,:) = cohen_cohere_results_high_low_mia(ii).baseline_subtracted_imaginary_coherence;
    low_high_coherence(ii,:,:) = cohen_cohere_results_low_high_mia(ii).imaginary_coherence;
    high_low_coherence(ii,:,:) = cohen_cohere_results_high_low_mia(ii).imaginary_coherence;
end

bsub_low_high_coherence = squeeze(mean(bs_low_high_coherence,1));
bsub_high_low_coherence = squeeze(mean(bs_high_low_coherence,1));

lh_imag_cohere = squeeze(mean(low_high_coherence));
hl_imag_cohere = squeeze(mean(high_low_coherence));


rdylbu = brewermap(1024,'-RdYlBu');
ylorrd = brewermap(1024,'YlOrRd');
figure("Position",[10 10 2000 600],"Visible","on")
subplot(1,2,1)
sgtitle('Baseline-subtracted coherence at transitions')
contourf(cohere_times,cohere_freqs,bsub_low_high_coherence,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('Low           Time from regime change (s)           High');
ylabel('Frequency (Hz)')
title('Low to High transitions')
colormap("jet")

subplot(1,2,2)
contourf(cohere_times,cohere_freqs,bsub_high_low_coherence,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('High           Time from regime change (s)           Low');
ylabel('Frequency (Hz)')
title('High to Low transitions')
c = colorbar();
ylabel(c,'Baseline subtracted coherence')
colormap("jet")

figure("Position",[10 10 2000 600],"Visible","on")
sgtitle('Coherence at transitions')
subplot(1,2,1)
contourf(cohere_times,cohere_freqs,lh_imag_cohere,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('Low           Time from regime change (s)           High');
ylabel('Frequency (Hz)')
colormap("jet")

subplot(1,2,2)
contourf(cohere_times,cohere_freqs,hl_imag_cohere,40,'linecolor','none')
set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
xlabel('High           Time from regime change (s)           Low');
ylabel('Frequency (Hz)')
c = colorbar();
ylabel(c,'Coherence');
colormap("jet")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PHASE LAG INDEX %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:size(cohen_pli_results_high_low_CTL,2)
    wpli_low_high_ctl(i,:,:) = cohen_pli_results_low_high_CTL(i).weighted_phase_lag_index_trials;
    wpli_high_low_ctl(i,:,:) = cohen_pli_results_high_low_CTL(i).weighted_phase_lag_index_trials;
    ispc_low_high_ctl(i,:,:) = cohen_pli_results_low_high_CTL(i).ispc_trials;
    ispc_high_low_ctl(i,:,:) = cohen_pli_results_high_low_CTL(i).ispc_trials;
end

wpli_hl_ctl = squeeze(nanmean(wpli_high_low_ctl,1));
wpli_lh_ctl = squeeze(nanmean(wpli_low_high_ctl,1));
ispc_hl_ctl = squeeze(nanmean(ispc_high_low_ctl,1));
ispc_lh_ctl = squeeze(nanmean(ispc_low_high_ctl,1));


freqs2use = cohen_pli_results_low_high_CTL(1).frequencies;
times2save = linspace(-1000,1000,2001);
ispc_times= cohen_pli_results_low_high_CTL(1).time;

figure("Position",[10 10 2000 600],"Visible","on")
subplot(1,2,1)
sgtitle('Weighted Phase-lag Index, Baseline Subtracted')
contourf(times2save(500:1501),freqs2use,imgaussfilt(wpli_lh_ctl(:,500:1501),1.2),'linecolor','none')
set(gca,'clim',[-.08 .08],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xlabel('Time from transition (msec)');
ylabel('Frequency (Hz)');
xline(0,'k--');
title('low to high transitions')
colormap("jet");
subplot(1,2,2)
contourf(times2save(500:1501),freqs2use,imgaussfilt(wpli_hl_ctl(:,500:1501),1.2),'linecolor','none')
set(gca,'clim',[-.08 .08],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xline(0,'k--');
xlabel('Time from transition (msec)');
ylabel('Frequency (Hz)');
title('high to low transitions')
c = colorbar();
ylabel(c,'PLI');
colormap("jet");


figure("Position",[10 10 2000 600],"Visible","on")
subplot(1,2,1)
sgtitle('Inter-site phase clustering, Baseline Subtracted')
contourf(times2save(500:1501),freqs2use,imgaussfilt(ispc_lh_ctl(:,500:1501),4),'linecolor','none')
set(gca,'clim',[-.05 .05],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xlabel('Time from transition (msec)');
ylabel('Frequency (Hz)');
xline(0,'k--');
title('low to high transitions')
colormap("jet");
subplot(1,2,2)
contourf(times2save(500:1501),freqs2use,imgaussfilt(ispc_hl_ctl(:,500:1501),4),'linecolor','none')
set(gca,'clim',[-.05 .05],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xline(0,'k--');
xlabel('Time from transition (msec)');
ylabel('Frequency (Hz)');
title('high to low transitions')
c = colorbar();
ylabel(c,'PLI');
colormap("jet");



for i = 1:size(cohen_pli_results_high_low_mia,2)
    wpli_low_high_mia(i,:,:) = cohen_pli_results_low_high_mia(i).weighted_phase_lag_index_trials;
    wpli_high_low_mia(i,:,:) = cohen_pli_results_high_low_mia(i).weighted_phase_lag_index_trials;
    ispc_low_high_mia(i,:,:) = cohen_pli_results_low_high_mia(i).ispc_trials;
    ispc_high_low_mia(i,:,:) = cohen_pli_results_high_low_mia(i).ispc_trials;
end

wpli_hl_mia = squeeze(nanmean(wpli_high_low_mia,1));
wpli_lh_mia = squeeze(nanmean(wpli_low_high_mia,1));
ispc_hl_mia = squeeze(nanmean(ispc_high_low_mia,1));
ispc_lh_mia = squeeze(nanmean(ispc_low_high_mia,1));

figure("Position",[10 10 2000 600],"Visible","on")
subplot(1,2,1)
sgtitle('Weighted Phase-lag Index, Baseline Subtracted')
contourf(times2save(500:1501),freqs2use,imgaussfilt(wpli_lh_mia(:,500:1501),1.2),'linecolor','none')
set(gca,'clim',[-.08 .08],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xlabel('Time from transition (msec)');
ylabel('Frequency (Hz)');
xline(0,'k--');
title('low to high transitions')
colormap("jet");
subplot(1,2,2)
contourf(times2save(500:1501),freqs2use,imgaussfilt(wpli_hl_mia(:,500:1501),1.2),'linecolor','none')
set(gca,'clim',[-.08 .08],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xline(0,'k--');
xlabel('Time from transition (msec)');
ylabel('Frequency (Hz)');
title('high to low transitions')
c = colorbar();
ylabel(c,'PLI');
colormap("jet");

figure("Position",[10 10 2000 600],"Visible","on")
subplot(1,2,1)
sgtitle('Inter-site phase clustering, Baseline Subtracted')
contourf(times2save(500:1501),freqs2use,imgaussfilt(ispc_lh_mia(:,500:1501),4),'linecolor','none')
set(gca,'clim',[-.05 .05],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xlabel('Time from transition (msec)');
ylabel('Frequency (Hz)');
xline(0,'k--');
title('low to high transitions')
colormap("jet");
subplot(1,2,2)
contourf(times2save(500:1501),freqs2use,imgaussfilt(ispc_hl_mia(:,500:1501),4),'linecolor','none')
set(gca,'clim',[-.05 .05],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xline(0,'k--');
xlabel('Time from transition (msec)');
ylabel('Frequency (Hz)');
title('high to low transitions')
c = colorbar();
ylabel(c,'PLI');
colormap("jet");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CROSS FREQUENCY COUPLING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('/Volumes/1TB ALPHA/MIA_Coherence_Results/CTL_master_cross_freq_results_v3.mat')
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/MIA_master_cross_freq_results_v3.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for z = 1:size(cross_freq_results_low_high_CTL,2)
    lh_ch_one_crf(z,:,:) = squeeze(mean(cross_freq_results_low_high_CTL(z).ch_one_cross_freq_coupling));
    lh_ch_two_crf(z,:,:) = squeeze(mean(cross_freq_results_low_high_CTL(z).ch_two_cross_freq_coupling));
    lh_ch_one_on_two_crf(z,:,:) = squeeze(mean(cross_freq_results_low_high_CTL(z).ch_one_on_two_cross_freq_coupling));
    lh_ch_two_one_one_crf(z,:,:) = squeeze(mean(cross_freq_results_low_high_CTL(z).ch_two_on_one_cross_freq_coupling));
end

freqVec1 = cross_freq_results_low_high_CTL.freqVec1;
freqVec2 = cross_freq_results_low_high_CTL.freqVec2;

ctl_lh_one(1:100,1:128) = squeeze(mean(lh_ch_one_crf,1));
ctl_lh_two(1:100,1:128) = squeeze(mean(lh_ch_two_crf,1));
ctl_lh_one_on_two(1:100,1:128) = squeeze(mean(lh_ch_one_on_two_crf,1));

contourf(freqVec1,freqVec2,ctl_lh_one,'linecolor','none')
contourf(freqVec1,freqVec2,ctl_lh_two,'linecolor','none')
contourf(freqVec1,freqVec2,ctl_lh_one_on_two,'linecolor','none')