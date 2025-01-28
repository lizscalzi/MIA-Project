
%%% Adapted from
%% Analyzing Neural Time Series Data
% Chapter 20, Mike X Cohen
% By Rob Munn, University of Otago, Department of Anatomy
% July 2024

%%% This code takes the struct EEG with fields
%%% EEG.data = a channels x time x trials matrix of data
%%% EEG.trials, EEG.pnts = arrays of number of trials, and number of data
%%% points
%%% EEG.srate = sample rate of the eeg
%%% also takes a descriptive field "animal name"
%%% returns a "results" struct with power, ERP power etc for each channel

function [results] = calculate_ERP_phase_locking(EEG,animal_name,plotfig)

% wavelet parameters
min_freq = 2;
max_freq = 60;
num_frex = 20;
%%% channels x points x trials
% baseline time window
baseline_time = [ -8000 -6000 ];
numb_msecs = (EEG.pnts-1)/EEG.srate * 1000;
EEG.times = linspace(-numb_msecs,numb_msecs,EEG.pnts);
% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution(1:2) = n_wavelet+n_data-1;
n_convolution(3)   = n_wavelet+EEG.pnts-1; % ERP is only one trial-length

for z = 1:2
    % find sensor index
    sensoridx = z;

    % compute ERP
    erp = squeeze(mean(EEG.data(sensoridx,:,:),3));

    % compute induced power by subtracting ERP from each trial
    induced_EEG = squeeze(EEG.data(sensoridx,:,:)) - repmat(erp',1,EEG.trials);

    % FFT of data
    fft_EEG{1} = fft(reshape(EEG.data(sensoridx,:,:),1,EEG.pnts*EEG.trials),n_convolution(1)); % total
    fft_EEG{2} = fft(reshape(induced_EEG,1,EEG.pnts*EEG.trials),n_convolution(2)); % induced
    fft_EEG{3} = fft(erp,n_convolution(3)); % evoked; note that it doesn't matter that the FFT is longer than the time series

    % convert baseline from ms to indices
    [~,baseidx(1)] = min(abs(EEG.times-baseline_time(1)));
    [~,baseidx(2)] = min(abs(EEG.times-baseline_time(2)));

    % initialize output time-frequency data
    tf = zeros(4,length(frequencies),EEG.pnts);

    for fi=1:length(frequencies)

        % create wavelet
        wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*(4/(2*pi*frequencies(fi)))^2))/frequencies(fi);

        % run convolution for each of total, induced, and evoked
        for i=1:3

            % take FFT of data
            fft_wavelet = fft(wavelet,n_convolution(i));

            % convolution...
            convolution_result_fft = ifft(fft_wavelet.*fft_EEG{i},n_convolution(i));
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);

            % reshaping and trial averaging is done only on all trials
            if i<3
                convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);

                % compute power
                tf(i,fi,:) = mean(abs(convolution_result_fft).^2,2);
            else
                % with only one trial-length, just compute power with no averaging
                tf(i,fi,:) = abs(convolution_result_fft).^2;
            end

            % db correct power
            tf(i,fi,:) = 10*log10( squeeze(tf(i,fi,:)) ./ mean(tf(i,fi,baseidx(1):baseidx(2)),3) );

            % inter-trial phase consistency on total EEG
            if i==1
                tf(4,fi,:) = abs(mean(exp(1i*angle(convolution_result_fft)),2));
            end
        end % end loop around total, evoked, induced
    end % end frequency loop

    % scale ERP for plotting
    erpt = (erp-min(erp))./max(erp-min(erp));
    erpt = erpt*(frequencies(end)-frequencies(1))+frequencies(1);

    if z == 1
        results.total_power_ch_one = squeeze(tf(1,:,:));
        results.non_phase_locked_power_ch_one = squeeze(tf(2,:,:));
        results.ERP_power_ch_one = squeeze(tf(3,:,:));
        results.ITPC_ch_one = squeeze(tf(4,:,:));
        results.ERP_ch_one = erpt;
    elseif z == 2
        results.total_power_ch_two = squeeze(tf(1,:,:));
        results.non_phase_locked_power_ch_two = squeeze(tf(2,:,:));
        results.ERP_power_ch_two = squeeze(tf(3,:,:));
        results.ITPC_ch_two = squeeze(tf(4,:,:));
        results.erp_ch_two = erpt;
    end
end

%%% CONTROLS
clear erp induced_EEG fft_EEG tf convolution_result_fft n_convolution n_data n_wavelet
EEG.trials = size(EEG.control_data,3);
% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution(1:2) = n_wavelet+n_data-1;
n_convolution(3)   = n_wavelet+EEG.pnts-1; % ERP is only one trial-length

for z = 1:2
    % find sensor index
    sensoridx = z;

    % compute ERP
    erp = squeeze(mean(EEG.control_data(sensoridx,:,:),3));

    % compute induced power by subtracting ERP from each trial
    induced_EEG = squeeze(EEG.control_data(sensoridx,:,:)) - repmat(erp',1,EEG.trials);

    % FFT of data
    fft_EEG{1} = fft(reshape(EEG.control_data(sensoridx,:,:),1,EEG.pnts*EEG.trials),n_convolution(1)); % total
    fft_EEG{2} = fft(reshape(induced_EEG,1,EEG.pnts*EEG.trials),n_convolution(2)); % induced
    fft_EEG{3} = fft(erp,n_convolution(3)); % evoked; note that it doesn't matter that the FFT is longer than the time series

    % convert baseline from ms to indices
    [~,baseidx(1)] = min(abs(EEG.times-baseline_time(1)));
    [~,baseidx(2)] = min(abs(EEG.times-baseline_time(2)));

    % initialize output time-frequency data
    tf = zeros(4,length(frequencies),EEG.pnts);

    for fi=1:length(frequencies)

        % create wavelet
        wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*(4/(2*pi*frequencies(fi)))^2))/frequencies(fi);

        % run convolution for each of total, induced, and evoked
        for i=1:3

            % take FFT of data
            fft_wavelet = fft(wavelet,n_convolution(i));

            % convolution...
            convolution_result_fft = ifft(fft_wavelet.*fft_EEG{i},n_convolution(i));
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);

            % reshaping and trial averaging is done only on all trials
            if i<3
                try
                convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
                catch
                    keyboard
                end

                % compute power
                tf(i,fi,:) = mean(abs(convolution_result_fft).^2,2);
            else
                % with only one trial-length, just compute power with no averaging
                tf(i,fi,:) = abs(convolution_result_fft).^2;
            end

            % db correct power
            tf(i,fi,:) = 10*log10( squeeze(tf(i,fi,:)) ./ mean(tf(i,fi,baseidx(1):baseidx(2)),3) );

            % inter-trial phase consistency on total EEG
            if i==1
                tf(4,fi,:) = abs(mean(exp(1i*angle(convolution_result_fft)),2));
            end
        end % end loop around total, evoked, induced
    end % end frequency loop

    % scale ERP for plotting
    erpt = (erp-min(erp))./max(erp-min(erp));
    erpt = erpt*(frequencies(end)-frequencies(1))+frequencies(1);
    if z == 1
        results.total_power_ch_one_control = squeeze(tf(1,:,:));
        results.non_phase_locked_power_ch_one_control = squeeze(tf(2,:,:));
        results.ERP_power_ch_one_control = squeeze(tf(3,:,:));
        results.ITPC_ch_one_control = squeeze(tf(4,:,:));
        results.ERP_ch_one_control = erpt;
    elseif z == 2
        results.total_power_ch_two_control = squeeze(tf(1,:,:));
        results.non_phase_locked_power_ch_two_control = squeeze(tf(2,:,:));
        results.ERP_power_ch_two_control = squeeze(tf(3,:,:));
        results.ITPC_ch_two_control = squeeze(tf(4,:,:));
        results.erp_ch_two_control = erpt;
    end
end



analysis_labels = {'Total';'Non-phase-locked';'ERP power';'ITPC'};
%
% % color limits
% clims = [ -3 3; -3 3; -12 12; 0 .6 ];
%
% % scale ERP for plotting
% erpt = (erp-min(erp))./max(erp-min(erp));
% erpt = erpt*(frequencies(end)-frequencies(1))+frequencies(1);
% plot_title = {'Channel one';'Channel two'};
%
% if plotfig == 1
%     figure(z)
%     for i=1:4
%         subplot(2,3,i)
%         contourf(EEG.times,frequencies,squeeze(tf(i,:,:)),40,'linecolor','none')
%         set(gca,'clim',clims(i,:),'xlim',[-1000 1000],'xtick',-1000:1000:800)
%         xlabel('Time (ms)')
%         ylabel('Frequency (Hz)')
%         title(analysis_labels{i})
%         % plot ERP on top
%         hold on
%         plot(EEG.times,erpt,'k')
%         sgtitle(plot_title{z})
%     end
%     savename = strcat(animal_name,'_',plot_title{z});
%     saveas(gcf,savename,'svg');
%     close all
% else
% end



results.freqs = frequencies;
results.times = EEG.times;
results.tf_types = analysis_labels;



end
