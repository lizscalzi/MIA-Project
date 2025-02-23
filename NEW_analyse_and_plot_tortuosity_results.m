%%%%%%%%%%%%%%%%%%%%
%% COHEN ANALYSES %%
%%%%%%%%%%%%%%%%%%%%

load('/Volumes/1TB ALPHA/MIA_Coherence_Results/CTL_master_cohen_results_v3.mat');
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/MIA_master_cohen_results_v3.mat');
plotfig = false;
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

if plotfig == true
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
else
end

index_to_one_second_around_transitions = 876:1126;
index_to_theta = 1:11;


%%% CHANNEL TWO
erp_data_channel_two_low_high_transitions(1,:,:) = squeeze(nanmean(total_low_high_channel_two));
erp_data_channel_two_low_high_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_low_high_channel_two));
erp_data_channel_two_low_high_transitions(3,:,:) = squeeze(nanmean(erp_power_low_high_channel_two));
erp_data_channel_two_low_high_transitions(4,:,:) = squeeze(nanmean(itpc_low_high_channel_two));

erp_control_two(1,:,:) = squeeze(nanmean(total_channel_two_control));
erp_control_two(2,:,:) = squeeze(nanmean(non_phase_locked_channel_two_control));
erp_control_two(3,:,:) = squeeze(nanmean(erp_power_channel_two_control));
erp_control_two(4,:,:) = squeeze(nanmean(itpc_channel_two_control));

if plotfig == true
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
else
end

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

if plotfig == true
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
        set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
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
        set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
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
else
end

erp_theta_range = 6:11;
for z = 1:size(total_low_high_channel_one,1)
    total_low_high_ch_one_theta_prechange(z,1) = mean(mean(total_low_high_channel_one(z,erp_theta_range,1:1000)));
    total_low_high_ch_one_theta_postchange(z,1) = mean(mean(total_low_high_channel_one(z,erp_theta_range,1001:2001)));
end

%%
[control_animal_power_results] = power_analysis_for_tortuosity(total_low_high_channel_one,total_low_high_channel_two,total_high_low_channel_one,total_high_low_channel_two,frequencies,times);
[control_animal_control_segments_power_results] = power_analysis_for_tortuosity(total_channel_one_control,total_channel_two_control,total_channel_one_control,total_channel_two_control,frequencies,times);
%%

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

if plotfig == true
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
else
end

%%% CHANNEL TWO
erp_data_channel_two_low_high_transitions(1,:,:) = squeeze(nanmean(total_low_high_channel_two));
erp_data_channel_two_low_high_transitions(2,:,:) = squeeze(nanmean(non_phase_locked_low_high_channel_two));
erp_data_channel_two_low_high_transitions(3,:,:) = squeeze(nanmean(erp_power_low_high_channel_two));
erp_data_channel_two_low_high_transitions(4,:,:) = squeeze(nanmean(itpc_low_high_channel_two));

erp_control_two(1,:,:) = squeeze(nanmean(total_channel_two_control));
erp_control_two(2,:,:) = squeeze(nanmean(non_phase_locked_channel_two_control));
erp_control_two(3,:,:) = squeeze(nanmean(erp_power_channel_two_control));
erp_control_two(4,:,:) = squeeze(nanmean(itpc_channel_two_control));

if plotfig == true
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
else
end

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

if plotfig == true
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
        set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
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
        set(gca,'clim',clims(i,:),'xlim',[-4 4],'xtick',-4:0.5:4)
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
else
end

erp_theta_range = 6:11;
for z = 1:size(total_low_high_channel_one,1)
    total_low_high_ch_one_theta_prechange(z,1) = mean(mean(total_low_high_channel_one(z,erp_theta_range,1:1000)));
    total_low_high_ch_one_theta_postchange(z,1) = mean(mean(total_low_high_channel_one(z,erp_theta_range,1001:2001)));
end


[mia_animal_power_results] = power_analysis_for_tortuosity(total_low_high_channel_one,total_low_high_channel_two,total_high_low_channel_one,total_high_low_channel_two,frequencies,times);
[mia_animal_control_segments_power_results] = power_analysis_for_tortuosity(total_channel_one_control,total_channel_two_control,total_channel_one_control,total_channel_two_control,frequencies,times);


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

if plotfig == true

    %%% raw imaginary coherence
    figure("Position",[10 10 2000 600],"Visible","on")
    sgtitle('Coherence at transitions')
    subplot(1,2,1)
    contourf(cohere_times,cohere_freqs,lh_imag_cohere,40,'linecolor','none')
    set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
    xlabel('Low           Time from regime change (s)           High');
    ylabel('Frequency (Hz)')
    xline(0,'k--');
    colormap("jet")
    subplot(1,2,2)
    contourf(cohere_times,cohere_freqs,hl_imag_cohere,40,'linecolor','none')
    set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
    xlabel('High           Time from regime change (s)           Low');
    ylabel('Frequency (Hz)')
    c = colorbar();
    ylabel(c,'Coherence');
    xline(0,'k--');
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
else
end


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

if plotfig == true
    figure()
    g = gramm('x',type,'y',low_high_control_theta_coherence,'group',axis);
    g.geom_line();
    g.set_order_options("x",0)
    g.update('x',type,'y',low_high_control_theta_coherence,'color',type)
    g.geom_point();
    g.set_names('x','','y','Baseline-subtracted coherence');
    g.set_title('Low to high transitions');
    g.axe_property('YLim',[0 0.1]);
    g.draw();

    figure()
    g = gramm('x',type,'y',high_low_control_theta_coherence,'group',axis);
    g.geom_line();
    g.set_order_options("x",0)
    g.update('x',type,'y',high_low_control_theta_coherence,'color',type)
    g.geom_point();
    g.set_names('x','','y','Baseline-subtracted coherence');
    g.set_title('Low to high transitions');
    g.axe_property('YLim',[-0.1 0.1]);
    g.draw();
else
end

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

if plotfig == true
    rdylbu = brewermap(1024,'-RdYlBu');
    ylorrd = brewermap(1024,'YlOrRd');
    figure("Position",[10 10 2000 600],"Visible","on")
    subplot(1,2,1)
    sgtitle('Baseline-subtracted coherence at transitions')
    contourf(cohere_times,cohere_freqs,bsub_low_high_coherence,40,'linecolor','none')
    set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
    xlabel('Low           Time from regime change (s)           High');
    xline(0,'k--');
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
    xline(0,'k--');
    ylabel(c,'Baseline subtracted coherence')
    colormap("jet")

    figure("Position",[10 10 2000 600],"Visible","on")
    sgtitle('Coherence at transitions')
    subplot(1,2,1)
    contourf(cohere_times,cohere_freqs,lh_imag_cohere,40,'linecolor','none')
    set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
    xlabel('Low           Time from regime change (s)           High');
    xline(0,'k--');
    ylabel('Frequency (Hz)')
    colormap("jet")

    subplot(1,2,2)
    contourf(cohere_times,cohere_freqs,hl_imag_cohere,40,'linecolor','none')
    set(gca,'clim',[0 0.17],'xlim',[-2 2],'ylim',[0 100],'xtick',-2:0.5:2,'yscale','log')
    xlabel('High           Time from regime change (s)           Low');
    ylabel('Frequency (Hz)')
    xline(0,'k--');
    c = colorbar();
    ylabel(c,'Coherence');
    colormap("jet")
else
end
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

if plotfig == true
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
else
end

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

if plotfig == true
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
else
end
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
    lh_ch_two_on_one_crf(z,:,:) = squeeze(mean(cross_freq_results_low_high_CTL(z).ch_two_on_one_cross_freq_coupling));
end

freqVec1 = cross_freq_results_low_high_CTL.freqVec1;
freqVec2 = cross_freq_results_low_high_CTL.freqVec2;

ctl_lh_one(1:100,1:128) = squeeze(mean(lh_ch_one_crf,1));
ctl_lh_two(1:100,1:128) = squeeze(mean(lh_ch_two_crf,1));
ctl_lh_one_on_two(1:100,1:128) = squeeze(mean(lh_ch_one_on_two_crf,1));
ctl_lh_two_on_one(1:100,1:128) = squeeze(mean(lh_ch_two_on_one_crf,1));

if plotfig == true
    figure()
    subplot(1,4,1)
    contourf(freqVec1,freqVec2,ctl_lh_one,'linecolor','none')
    subplot(1,4,2)
    contourf(freqVec1,freqVec2,ctl_lh_two,'linecolor','none')
    subplot(1,4,3)
    contourf(freqVec1,freqVec2,ctl_lh_one_on_two,'linecolor','none')
    subplot(1,4,4)
    contourf(freqVec1,freqVec2,ctl_lh_two_on_one,'linecolor','none')
    colormap('jet')
else
end

%% BY SESSION POWER RESULTS AND PLOTTING

[freq_power_stats] = by_session_power_results_and_plotting(control_animal_power_results,control_animal_control_segments_power_results,mia_animal_power_results,mia_animal_control_segments_power_results,plotfig,frequencies,times);


%% PHASE

load('/Volumes/1TB ALPHA/MIA_Coherence_Results/CTL_master_phase_results_v3.mat');
load('/Volumes/1TB ALPHA/MIA_Coherence_Results/MIA_master_phase_results_v3.mat');


%%% CONTROL
for z = 1:size(phase_results_CTL,2)
    psi_theta_lh(:,z) = nanmean(phase_results_CTL(z).low_high.psi_theta_hpc_to_pfc,2);
    psi_theta_hl(:,z) = nanmean(phase_results_CTL(z).high_low.psi_theta_hpc_to_pfc,2);
    psi_theta_control(:,z) = nanmean(phase_results_CTL(z).control.psi_theta_hpc_to_pfc,2);

    psi_slow_gamma_lh(:,z) = nanmean(phase_results_CTL(z).low_high.psi_slow_gamma_hpc_to_pfc,2);
    psi_slow_gamma_hl(:,z) = nanmean(phase_results_CTL(z).high_low.psi_slow_gamma_hpc_to_pfc,2);
    psi_slow_gamma_control(:,z) = nanmean(phase_results_CTL(z).control.psi_slow_gamma_hpc_to_pfc,2);

    psi_fast_gamma_lh(:,z) = nanmean(phase_results_CTL(z).low_high.psi_fast_gamma_hpc_to_pfc,2);
    psi_fast_gamma_hl(:,z) = nanmean(phase_results_CTL(z).high_low.psi_fast_gamma_hpc_to_pfc,2);
    psi_fast_gamma_control(:,z) = nanmean(phase_results_CTL(z).control.psi_fast_gamma_hpc_to_pfc,2);
end
window = [1 125 250 375 500 625 750 875 1000 1125 1250 1375 1500 1625 1750 1875 2000];
window = (window/250)-4;

figure()
shadedErrorBar(window(1:16),nanmean(psi_theta_lh,2),nanstd(psi_theta_lh,0,2)./sqrt(size(psi_theta_lh,2)),'lineProps','k--');
hold on
shadedErrorBar(window(1:16),nanmean(psi_slow_gamma_lh,2),nanstd(psi_slow_gamma_lh,0,2)./sqrt(size(psi_slow_gamma_lh,2)),'lineProps','b-');
shadedErrorBar(window(1:16),nanmean(psi_fast_gamma_lh,2),nanstd(psi_fast_gamma_lh,0,2)./sqrt(size(psi_fast_gamma_lh,2)),'lineProps','r-');
xlabel('Time from switch (s)');
ylabel('Phase slope index');
legend('Theta (5:12Hz)','[]','[]','[]','Slow Gamma (25-50Hz)','[]','[]','[]','Fast Gamma (65-125Hz)')
yline(0,'k--');
xline(0,'k--');
title('PSI, low to high transitions')
hold off

figure()
shadedErrorBar(window(1:16),nanmean(psi_theta_hl,2),nanstd(psi_theta_hl,0,2)./sqrt(size(psi_theta_hl,2)),'lineProps','k--');
hold on
shadedErrorBar(window(1:16),nanmean(psi_slow_gamma_hl,2),nanstd(psi_slow_gamma_hl,0,2)./sqrt(size(psi_slow_gamma_hl,2)),'lineProps','b-');
shadedErrorBar(window(1:16),nanmean(psi_fast_gamma_hl,2),nanstd(psi_fast_gamma_hl,0,2)./sqrt(size(psi_fast_gamma_hl,2)),'lineProps','r-');
xlabel('Time from switch (s)');
ylabel('Phase slope index');
legend('Theta (5:12Hz)','[]','[]','[]','Slow Gamma (25-50Hz)','[]','[]','[]','Fast Gamma (65-125Hz)')
yline(0,'k--');
xline(0,'k--');
title('PSI, high to low transitions');
hold off


%%% MIA
for z = 1:size(phase_results_MIA,2)
    psi_theta_lh_mia(:,z) = nanmean(phase_results_MIA(z).low_high.psi_theta_hpc_to_pfc,2);
    psi_theta_hl_mia(:,z) = nanmean(phase_results_MIA(z).high_low.psi_theta_hpc_to_pfc,2);
    psi_theta_control_mia(:,z) = nanmean(phase_results_MIA(z).control.psi_theta_hpc_to_pfc,2);

    psi_slow_gamma_lh_mia(:,z) = nanmean(phase_results_MIA(z).low_high.psi_slow_gamma_hpc_to_pfc,2);
    psi_slow_gamma_hl_mia(:,z) = nanmean(phase_results_MIA(z).high_low.psi_slow_gamma_hpc_to_pfc,2);
    psi_slow_gamma_control_mia(:,z) = nanmean(phase_results_MIA(z).control.psi_slow_gamma_hpc_to_pfc,2);

    psi_fast_gamma_lh_mia(:,z) = nanmean(phase_results_MIA(z).low_high.psi_fast_gamma_hpc_to_pfc,2);
    psi_fast_gamma_hl_mia(:,z) = nanmean(phase_results_MIA(z).high_low.psi_fast_gamma_hpc_to_pfc,2);
    psi_fast_gamma_control_mia(:,z) = nanmean(phase_results_MIA(z).control.psi_fast_gamma_hpc_to_pfc,2);
end
window = [1 125 250 375 500 625 750 875 1000 1125 1250 1375 1500 1625 1750 1875 2000];
window = (window/250)-4;

figure()
shadedErrorBar(window(1:16),nanmean(psi_theta_lh_mia,2),nanstd(psi_theta_lh_mia,0,2)./sqrt(size(psi_theta_lh_mia,2)),'lineProps','k--');
hold on
shadedErrorBar(window(1:16),nanmean(psi_slow_gamma_lh_mia,2),nanstd(psi_slow_gamma_lh_mia,0,2)./sqrt(size(psi_slow_gamma_lh_mia,2)),'lineProps','b-');
shadedErrorBar(window(1:16),nanmean(psi_fast_gamma_lh_mia,2),nanstd(psi_fast_gamma_lh_mia,0,2)./sqrt(size(psi_fast_gamma_lh_mia,2)),'lineProps','r-');
xlabel('Time from switch (s)');
ylabel('Phase slope index');
legend('Theta (5:12Hz)','[]','[]','[]','Slow Gamma (25-50Hz)','[]','[]','[]','Fast Gamma (65-125Hz)')
yline(0,'k--');
xline(0,'k--');
title('PSI, low to high transitions; MIA animals')
hold off

figure()
shadedErrorBar(window(1:16),nanmean(psi_theta_hl_mia,2),nanstd(psi_theta_hl_mia,0,2)./sqrt(size(psi_theta_hl_mia,2)),'lineProps','k--');
hold on
shadedErrorBar(window(1:16),nanmean(psi_slow_gamma_hl_mia,2),nanstd(psi_slow_gamma_hl_mia,0,2)./sqrt(size(psi_slow_gamma_hl_mia,2)),'lineProps','b-');
shadedErrorBar(window(1:16),nanmean(psi_fast_gamma_hl_mia,2),nanstd(psi_fast_gamma_hl_mia,0,2)./sqrt(size(psi_fast_gamma_hl_mia,2)),'lineProps','r-');
xlabel('Time from switch (s)');
ylabel('Phase slope index');
legend('Theta (5:12Hz)','[]','[]','[]','Slow Gamma (25-50Hz)','[]','[]','[]','Fast Gamma (65-125Hz)')
yline(0,'k--');
xline(0,'k--');
title('PSI, high to low transitions; MIA animals');
hold off
