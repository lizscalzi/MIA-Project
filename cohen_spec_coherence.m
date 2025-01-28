function [results] = cohen_spec_coherence(EEG,animal_name,plotfig)


freqs2use  = logspace(log10(2),log10(120),118); % 2-120 Hz in 118 steps
times2save = -4000:50:4000;
% wavelet and FFT parameters
numb_msecs = (EEG.pnts-1)/EEG.srate * 1000;
EEG.times = linspace(-numb_msecs,numb_msecs,EEG.pnts);
time          = -1:1/EEG.srate:1;
half_wavelet  = (length(time)-1)/2;
num_cycles    = logspace(log10(4),log10(8),length(freqs2use));
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
baselinetm = [-7000 -6000];

chanidx    = zeros(1,2); % always initialize!
chanidx(1) = 1;
chanidx(2) = 2;

% data FFTs
data_fft1 = fft(reshape(EEG.data(chanidx(1),:,:),1,n_data),n_convolution);
data_fft2 = fft(reshape(EEG.data(chanidx(2),:,:),1,n_data),n_convolution);

times2saveidx = dsearchn(EEG.times',times2save');
baselineidx   = dsearchn(times2save',baselinetm');
% initialize
spectcoher = zeros(length(freqs2use),length(times2save));

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

    % compute power and cross-spectral power
    spec1 = mean(sig1.*conj(sig1),2);
    spec2 = mean(sig2.*conj(sig2),2);
    specX = abs(mean(sig1.*conj(sig2),2)).^2;

    % alternative notation for the same procedure, using the Euler-like expression: Me^ik
    %spec1 = mean(abs(sig1).^2,2);
    %spec2 = mean(abs(sig2).^2,2);
    %specX = abs(mean( abs(sig1).*abs(sig2) .* exp(1i*(angle(sig1)-angle(sig2))) ,2)).^2;

    % compute spectral coherence, using only requested time points
    spectcoher(fi,:) = specX(times2saveidx)./(spec1(times2saveidx).*spec2(times2saveidx));

    % yet another equivalent notation, just FYI
    %spec1 = sum(sig1.*conj(sig1),2);
    %spec2 = sum(sig2.*conj(sig2),2);
    %specX = sum(sig1.*conj(sig2),2);
    %spectcoher(fi,:) = abs(specX(times2saveidx)./sqrt(spec1(times2saveidx).*spec2(times2saveidx))).^2;


    % imaginary coherence
    spec1 = sum(sig1.*conj(sig1),2);
    spec2 = sum(sig2.*conj(sig2),2);
    specX = sum(sig1.*conj(sig2),2);
    imagspectcoher(fi,:) = abs(imag(specX(times2saveidx)./sqrt(spec1(times2saveidx).*spec2(times2saveidx))));
end

baseline_subtracted_cohere = spectcoher-repmat(mean(spectcoher(:,baselineidx(1):baselineidx(2)),2),1,size(spectcoher,2));
baseline_subtracted_imag_cohere = imagspectcoher-repmat(mean(spectcoher(:,baselineidx(1):baselineidx(2)),2),1,size(spectcoher,2));

results.coherence = spectcoher;
results.baseline_subtracted_coherence = baseline_subtracted_cohere;
results.imaginary_coherence = imagspectcoher;
results.baseline_subtracted_imaginary_coherence = baseline_subtracted_imag_cohere;
results.freqs = freqs2use;
results.time = times2save;

if plotfig == 1
figure("Position",[10 10 2000 600],"Visible","off");
subplot(121)
contourf(times2save,freqs2use,spectcoher,20,'linecolor','none') %
set(gca,'clim',[0 round(max(spectcoher,[],'all'),1)],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)),'xlim',[times2save(1) times2save(end)])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
hcb = colorbar;
ylabel(hcb,"Coherence");
title('"Raw" spectral coherence')
xline(0,'k--','LineWidth',4);
subplot(122)
contourf(times2save,freqs2use,baseline_subtracted_cohere,20,'linecolor','none') %
set(gca,'clim',[round(min(baseline_subtracted_cohere,[],'all'),1) round(max(baseline_subtracted_cohere,[],'all'),1)],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)),'xlim',[times2save(1) times2save(end)])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
hcb = colorbar;
ylabel(hcb,"Delta Coherence");
xline(0,'k--','LineWidth',4);
colormap("jet")
title('Baseline-subtracted spectral coherence')
savename = strcat(animal_name,'_','coherence_plot');
saveas(gcf,savename,'svg');
close all


figure("Position",[10 10 2000 600],"Visible","off");
subplot(121)
contourf(times2save,freqs2use,imagspectcoher,20,'linecolor','none') %
set(gca,'clim',[0 round(max(imagspectcoher,[],'all'),1)],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)),'xlim',[times2save(1) times2save(end)])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
hcb = colorbar;
ylabel(hcb,"Coherence");
title('"Raw" imaginary spectral coherence')
xline(0,'k--','LineWidth',4);
subplot(122)
contourf(times2save,freqs2use,baseline_subtracted_imag_cohere,20,'linecolor','none') %
set(gca,'clim',[round(min(baseline_subtracted_imag_cohere,[],'all'),1) round(max(baseline_subtracted_imag_cohere,[],'all'),1)],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)),'xlim',[times2save(1) times2save(end)])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
hcb = colorbar;
ylabel(hcb,"Delta Coherence");
title('Baseline-subtracted imaginary spectral coherence')
xline(0,'k--','LineWidth',4);
colormap("jet")
savename = strcat(animal_name,'_','imaginary_coherence_plot');
saveas(gcf,savename,'svg');
close all
else
end