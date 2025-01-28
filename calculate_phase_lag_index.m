
%%% Adapted from
%% Analyzing Neural Time Series Data
%% Mike X Cohen
% By Rob Munn, University of Otago, Department of Anatomy
% July 2024

%%% This code takes the struct EEG with fields
%%% EEG.data = a channels x time x trials matrix of data
%%% EEG.trials, EEG.pnts = arrays of number of trials, and number of data
%%% points
%%% EEG.srate = sample rate of the eeg
%%% also takes a descriptive field "animal name"
%%% returns a "results" struct with the phase lag index and intertrial
%%% phase clustering over time and over trials



function [results] = calculate_phase_lag_index(EEG,animal_name,plotfig)

% select channels
channel1 = 'hpc';
channel2 = 'pfc';

% specify some time-frequency parameters
freqs2use  = logspace(log10(2),log10(120),25); % 2-120 Hz in 25 steps
times2save = -1000:10:1000;
timewindow = linspace(1.5,3,length(freqs2use)); % number of cycles on either end of the center point (1.5 means a total of 3 cycles))
baselinetm = [-8000 -7000];
numb_msecs = (EEG.pnts-1)/EEG.srate * 1000;
EEG.times = linspace(-numb_msecs,numb_msecs,EEG.pnts);
% wavelet and FFT parameters
time          = -1:1/EEG.srate:1;
half_wavelet  = (length(time)-1)/2;
num_cycles    = logspace(log10(4),log10(8),length(freqs2use));
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

% time in indices
times2saveidx = dsearchn(EEG.times',times2save');
baselineidxF  = dsearchn(EEG.times',baselinetm');  % for the full temporal resolution data (thanks to Daniel Roberts for finding/reporting this bug here!)
baselineidx   = dsearchn(times2save',baselinetm'); % for the temporally downsampled data

chanidx = zeros(1,2); % always initialize!
chanidx(1) = 1;
chanidx(2) = 2;

% data FFTs
data_fft1 = fft(reshape(EEG.data(chanidx(1),:,:),1,n_data),n_convolution);
data_fft2 = fft(reshape(EEG.data(chanidx(2),:,:),1,n_data),n_convolution);

% initialize
ispc    = zeros(length(freqs2use),EEG.pnts);
pli     = zeros(length(freqs2use),EEG.pnts);
wpli    = zeros(length(freqs2use),EEG.pnts);
dwpli   = zeros(length(freqs2use),EEG.pnts);
dwpli_t = zeros(length(freqs2use),length(times2save));
ispc_t  = zeros(length(freqs2use),length(times2save));

for fi=1:length(freqs2use)
    
    % create wavelet and take FFT
    s = num_cycles(fi)/(2*pi*freqs2use(fi));
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig1 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig2 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % cross-spectral density
    cdd = sig1 .* conj(sig2);
    
    % ISPC
    ispc(fi,:) = abs(mean(exp(1i*angle(cdd)),2)); % note: equivalent to ispc(fi,:) = abs(mean(exp(1i*(angle(sig1)-angle(sig2))),2));
    
    
    % take imaginary part of signal only
    cdi = imag(cdd);
    
    % phase-lag index
    pli(fi,:)  = abs(mean(sign(imag(cdd)),2));
    
    % weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
    wpli(fi,:) = abs( mean( abs(cdi).*sign(cdi) ,2) )./mean(abs(cdi),2);
    
    % debiased weighted phase-lag index (shortcut, as implemented in fieldtrip)
    imagsum      = sum(cdi,2);
    imagsumW     = sum(abs(cdi),2);
    debiasfactor = sum(cdi.^2,2);
    dwpli(fi,:)  = (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor);
    
    % compute time window in indices for this frequency
    time_window_idx = round((1000/freqs2use(fi))*timewindow(fi)/(1000/EEG.srate));

    for ti=1:length(times2save)
        imagsum        = sum(cdi(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:),1);
        imagsumW       = sum(abs(cdi(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:)),1);
        debiasfactor   = sum(cdi(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:).^2,1);
        dwpli_t(fi,ti) = mean((imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor));

        % compute phase synchronization
        phasesynch     = abs(mean(exp(1i*angle(cdd(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:))),1));
        ispc_t(fi,ti)  = mean(phasesynch);
    end
end

% baseline subtraction from all measures
ispc    = bsxfun(@minus,ispc,mean(ispc(:,baselineidxF(1):baselineidxF(2)),2));
ispc_t  = bsxfun(@minus,ispc_t,mean(ispc_t(:,baselineidx(1):baselineidx(2)),2));
pli     = bsxfun(@minus,pli,mean(pli(:,baselineidxF(1):baselineidxF(2)),2));
dwpli   = bsxfun(@minus,dwpli,mean(dwpli(:,baselineidxF(1):baselineidxF(2)),2));
dwpli_t = bsxfun(@minus,dwpli_t,mean(dwpli_t(:,baselineidx(1):baselineidx(2)),2));

if plotfig == 1
figure("Position",[10 10 2000 600],"Visible","off")
subplot(221)
contourf(times2save,freqs2use,pli(:,times2saveidx),40,'linecolor','none')
set(gca,'clim',[-.3 .3],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
title('PLI over trials')
colorbar()
subplot(222)
contourf(times2save,freqs2use,dwpli(:,times2saveidx),40,'linecolor','none')
set(gca,'clim',[-.2 .2],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
title('dWPLI over trials')
colorbar()
subplot(223)
contourf(times2save,freqs2use,ispc_t,40,'linecolor','none')
set(gca,'clim',[-.1 .1],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
title('ISPC over time')
colorbar()
subplot(224)
contourf(times2save,freqs2use,dwpli_t,40,'linecolor','none')
 set(gca,'clim',[-.1 .1],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
title('dWPLI over time')
colorbar()
colormap("jet");
animal_name_split = strsplit(animal_name,',');
cd('/Volumes/1TB ALPHA/MIA_Coherence_Results')
savename = strcat('/Volumes/1TB ALPHA/MIA_Coherence_Results/',animal_name_split{1},'_phase_lag_index');
saveas(gcf,savename,'svg');
else
end

results.ispc_trials = ispc;
results.phase_lag_index_trials = pli;
results.weighted_phase_lag_index_trials = dwpli;
results.ispc_time = ispc_t;
results.weighted_phase_lag_index_time = dwpli_t;
results.time = times2save;
results.time_index = times2saveidx;
results.frequencies = freqs2use;
results.animal = animal_name;

end