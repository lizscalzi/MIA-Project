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
%%% returns a "results" struct with the correlation of total power and
%%% spectral power between two brain regions centered on transitions or
%%% events

function [results] = calculate_power_correlation(EEG,animal_name,plotfig)

animal_name_split = strsplit(animal_name,',');
sensor1 = {'hippocampus'};
sensor2 = {'PFC'};
numb_msecs = (EEG.pnts-1)/EEG.srate * 1000;
EEG.times = linspace(-numb_msecs,numb_msecs,EEG.pnts);

timewin1 = [ -300 -100 ]; % in ms relative to stim onset
timewin2 = [  200  400 ];

centerfreq1 =  6; % in Hz
centerfreq2 =  6;


% convert time from ms to index
timeidx1 = zeros(size(timewin1));
timeidx2 = zeros(size(timewin2));
for i=1:2
    [~,timeidx1(i)] = min(abs(EEG.times-timewin1(i)));
    [~,timeidx2(i)] = min(abs(EEG.times-timewin2(i)));
end

% setup wavelet convolution and outputs
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
wavelet_cycles= 4.5;

% FFT of data (note: this doesn't change on frequency iteration)
fft_data1 = fft(reshape(EEG.data(1,:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_data2 = fft(reshape(EEG.data(2,:,:),1,EEG.pnts*EEG.trials),n_convolution);

% initialize output time-frequency data
%corrdata = zeros(EEG.trials,2);

% create wavelet and run convolution
fft_wavelet            = fft(exp(2*1i*pi*centerfreq1.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq1))^2)),n_convolution);
convolution_result_fft = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq1));
convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
analyticsignal1        = abs(convolution_result_fft).^2;

fft_wavelet            = fft(exp(2*1i*pi*centerfreq2.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq2))^2)),n_convolution);
convolution_result_fft = ifft(fft_wavelet.*fft_data2,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq2));
convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
analyticsignal2        = abs(convolution_result_fft).^2;


% Panel A: correlation in a specified window
tfwindowdata1 = mean(analyticsignal1(timeidx1(1):timeidx1(2),:),1);
tfwindowdata2 = mean(analyticsignal2(timeidx2(1):timeidx2(2),:),1);
% if plotfig == 1
% figure
% subplot(121)
% plot(tfwindowdata1,tfwindowdata2,'.')
% axis square
% title([ 'TF window correlation, r_p=' num2str(corr(tfwindowdata1',tfwindowdata2','type','p')) ])
% xlabel([ sensor1 ': ' num2str(timewin1(1)) '-' num2str(timewin1(2)) '; ' num2str(centerfreq1) ' Hz' ])
% ylabel([ sensor2 ': ' num2str(timewin2(1)) '-' num2str(timewin2(2)) '; ' num2str(centerfreq2) ' Hz' ])
% 
% 
% % also plot rank-transformed data
% subplot(122)
% plot(tiedrank(tfwindowdata1),tiedrank(tfwindowdata2),'.')
% axis square
% xlabel([ sensor1 ': ' num2str(timewin1(1)) '-' num2str(timewin1(2)) '; ' num2str(centerfreq1) ' Hz' ])
% ylabel([ sensor2 ': ' num2str(timewin2(1)) '-' num2str(timewin2(2)) '; ' num2str(centerfreq2) ' Hz' ])
% title([ 'TF window correlation, r_p=' num2str(corr(tfwindowdata1',tfwindowdata2','type','s')) ])
% else
% end

% panel B: correlation over time
corr_ts = zeros(size(EEG.times));
for ti=1:EEG.pnts
    corr_ts(ti) = corr(analyticsignal1(ti,:)',analyticsignal2(ti,:)','type','s');
end

if plotfig == 1
figure("Visible","off")
plot(EEG.times,corr_ts)
set(gca,'xlim',[-1000 1500])
xlabel('Time (ms)'), ylabel('Spearman''s rho')
savename = strcat('/Volumes/1TB ALPHA/MIA_Coherence_Results/',animal_name_split{1},'_power_correlation');
saveas(gcf,savename,'svg');
else 
end
results.power_correlation_time_over_trials = corr_ts(1,876:1188);
results.power_correlation_times = EEG.times(876:1188);

% Panel C: exploratory time-frequency power correlations
times2save = -1000:50:1500;
frex = logspace(log10(2),log10(120),40);


times2save_idx = zeros(size(times2save));
for i=1:length(times2save)
    [~,times2save_idx(i)] = min(abs(EEG.times-times2save(i)));
end

% rank-transforming the data can happen outside the frequency loop
seeddata_rank = tiedrank(tfwindowdata2);

% initialize output correlation matrix
expl_corrs = zeros(length(frex),length(times2save));

for fi=1:length(frex)
    
    % get power (via wavelet convolution) from signal1
    fft_wavelet            = fft(exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frex(fi)))^2)),n_convolution);
    convolution_result_fft = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*frex(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    analyticsignal1        = abs(convolution_result_fft).^2;
    
    for ti=1:length(times2save)
        expl_corrs(fi,ti) = 1-6*sum((seeddata_rank-tiedrank(analyticsignal1(times2save_idx(ti),:))).^2)/(EEG.trials*(EEG.trials^2-1));
    end
end
if plotfig == 1
figure("Position",[10 10 2000 600],"Visible","off")
contourf(times2save,frex,expl_corrs,40,'linecolor','none')
set(gca,'clim',[-.4 .4],'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),8)))
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Correlation over trials from seed ' sensor2 ', ' num2str(centerfreq2) ' Hz and ' num2str(timewin2(1)) '-' num2str(timewin2(2)) ' ms' ])
colorbar
colormap("jet")
savename = strcat('/Volumes/1TB ALPHA/MIA_Coherence_Results/',animal_name_split{1},'_spectral_power_correlation');
saveas(gcf,savename,'svg');
else 
end

results.spectral_power_correlation_over_trials = expl_corrs;
results.spectral_power_times = times2save;
results.spectral_power_frequencies = frex;


end