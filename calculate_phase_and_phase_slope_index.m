function [results] = calculate_phase_and_phase_slope_index(lh_eeg,hl_eeg)


%% PHASE SLOPE INDEX
freqbins(1,1) = 5; freqbins(1,2) = 12;
freqbins(2,1) = 25; freqbins(2,2) = 50;
freqbins(3,1) = 65; freqbins(3,2) = 125;


times = linspace(-4,4,2001);
results.times = times;


data(:,:,:) = lh_eeg.data(:,1:2000,:);
num_trials = size(data,3);

window = [1 125 250 375 500 625 750 875 1000 1125 1250 1375 1500 1625 1750 1875 2000];
results.window = window;
for h = 1:num_trials
for i = 1:size(window,2)-1
[out] = phase_slope_index(data(:,window(i)+1:window(i+1),h),250,freqbins,0);
results.low_high.psi_theta_hpc_to_pfc(i,h) = out(1,2,1);
results.low_high.psi_theta_pfc_to_hpc(i,h) = out(2,1,1);

results.low_high.psi_slow_gamma_hpc_to_pfc(i,h) = out(1,2,2);
results.low_high.psi_slow_gamma_pfc_to_hpc(i,h) = out(2,1,2);

results.low_high.psi_fast_gamma_hpc_to_pfc(i,h) = out(1,2,3);
results.low_high.psi_fast_gamma_pfc_to_hpc(i,h) = out(2,1,3);
clear out
end
end

clear data
data(:,:,:) = hl_eeg.data(:,1:2000,:);
num_trials = size(data,3);
for h = 1:num_trials
for i = 1:size(window,2)-1
[out] = phase_slope_index(data(:,window(i)+1:window(i+1),h),250,freqbins,0);
results.high_low.psi_theta_hpc_to_pfc(i,h) = out(1,2,1);
results.high_low.psi_theta_pfc_to_hpc(i,h) = out(2,1,1);

results.high_low.psi_slow_gamma_hpc_to_pfc(i,h) = out(1,2,2);
results.high_low.psi_slow_gamma_pfc_to_hpc(i,h) = out(2,1,2);

results.high_low.psi_fast_gamma_hpc_to_pfc(i,h) = out(1,2,3);
results.high_low.psi_fast_gamma_pfc_to_hpc(i,h) = out(2,1,3);
clear out
end
end


clear data
data(:,:,:) = lh_eeg.control_data(:,1:2000,:);
num_trials = size(data,3);
for h = 1:num_trials
for i = 1:size(window,2)-1
[out] = phase_slope_index(data(:,window(i)+1:window(i+1),h),250,freqbins,0);
results.control.psi_theta_hpc_to_pfc(i,h) = out(1,2,1);
results.control.psi_theta_pfc_to_hpc(i,h) = out(2,1,1);

results.control.psi_slow_gamma_hpc_to_pfc(i,h) = out(1,2,2);
results.control.psi_slow_gamma_pfc_to_hpc(i,h) = out(2,1,2);

results.control.psi_fast_gamma_hpc_to_pfc(i,h) = out(1,2,3);
results.control.psi_fast_gamma_pfc_to_hpc(i,h) = out(2,1,3);
clear out
end
end


%% PHASE OFFSET
%%% High to Low %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hpc_hl_eeg = squeeze(hl_eeg.data(1,:,:));
pfc_hl_eeg = squeeze(hl_eeg.data(2,:,:));

%%%% theta phase
for z = 1:size(hpc_hl_eeg,2)
    theta_hl_hpc = bandpass_colgin(hpc_hl_eeg(:,z),5,12,250);
    theta_hl_hpc = hilbert(theta_hl_hpc);
    hpc_phase_hl(z,:) = angle(theta_hl_hpc);
    theta_hl_pfc = bandpass_colgin(pfc_hl_eeg(:,z),5,12,250);
    theta_hl_pfc = hilbert(theta_hl_pfc);
    pfc_phase_hl(z,:) = angle(theta_hl_pfc);
    hl_theta_ang_diff(z,:) = angdiff(hpc_phase_hl(z,:),pfc_phase_hl(z,:));  
end

results.high_low.theta_phase_hpc = hpc_phase_hl;
results.high_low.theta_phase_pfc = pfc_phase_hl;
results.high_low.theta_angle_diff = hl_theta_ang_diff;
clear hpc_phase_hl pfc_phase_hl hl_theta_ang_diff

%%%% slow gamma phase
for z = 1:size(hpc_hl_eeg,2)
    slow_gamma_hl_hpc = bandpass_colgin(hpc_hl_eeg(:,z),25,50,250);
    slow_gamma_hl_hpc = hilbert(slow_gamma_hl_hpc);
    hpc_phase_hl(z,:) = angle(slow_gamma_hl_hpc);
    slow_gamma_hl_pfc = bandpass_colgin(pfc_hl_eeg(:,z),5,12,250);
    slow_gamma_hl_pfc = hilbert(slow_gamma_hl_pfc);
    pfc_phase_hl(z,:) = angle(slow_gamma_hl_pfc);
    hl_slow_gamma_ang_diff(z,:) = angdiff(hpc_phase_hl(z,:),pfc_phase_hl(z,:));
end

results.high_low.slow_gamma_phase_hpc = hpc_phase_hl;
results.high_low.slow_gamma_phase_pfc = pfc_phase_hl;
results.high_low.slow_gamma_angle_diff = hl_slow_gamma_ang_diff;
clear hpc_phase_hl pfc_phase_hl hl_slow_gamma_ang_diff

%%%% fast gamma phase
for z = 1:size(hpc_hl_eeg,2)
    fast_gamma_hl_hpc = bandpass_colgin(hpc_hl_eeg(:,z),65,125,250);
    fast_gamma_hl_hpc = hilbert(fast_gamma_hl_hpc);
    hpc_phase_hl(z,:) = angle(fast_gamma_hl_hpc);
    fast_gamma_hl_pfc = bandpass_colgin(pfc_hl_eeg(:,z),5,12,250);
    fast_gamma_hl_pfc = hilbert(fast_gamma_hl_pfc);
    pfc_phase_hl(z,:) = angle(fast_gamma_hl_pfc);
    hl_fast_gamma_ang_diff(z,:) = angdiff(hpc_phase_hl(z,:),pfc_phase_hl(z,:));
end

results.high_low.fast_gamma_phase_hpc = hpc_phase_hl;
results.high_low.fast_gamma_phase_pfc = pfc_phase_hl;
results.high_low.fast_gamma_angle_diff = hl_fast_gamma_ang_diff;
clear hpc_phase_hl pfc_phase_hl hl_fast_gamma_ang_diff


%%% Low to High %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hpc_lh_eeg = squeeze(lh_eeg.data(1,:,:));
pfc_lh_eeg = squeeze(lh_eeg.data(2,:,:));

%%%% theta phase
for z = 1:size(hpc_lh_eeg,2)
    theta_lh_hpc = bandpass_colgin(hpc_lh_eeg(:,z),5,12,250);
    theta_lh_hpc = hilbert(theta_lh_hpc);
    hpc_phase_lh(z,:) = angle(theta_lh_hpc);
    theta_lh_pfc = bandpass_colgin(pfc_lh_eeg(:,z),5,12,250);
    theta_lh_pfc = hilbert(theta_lh_pfc);
    pfc_phase_lh(z,:) = angle(theta_lh_pfc);
    lh_theta_ang_diff(z,:) = angdiff(hpc_phase_lh(z,:),pfc_phase_lh(z,:));
end

results.low_high.theta_phase_hpc = hpc_phase_lh;
results.low_high.theta_phase_pfc = pfc_phase_lh;
results.low_high.theta_angle_diff = lh_theta_ang_diff;
clear hpc_phase_lh pfc_phase_lh lh_theta_ang_diff



%%%% slow gamma phase
for z = 1:size(hpc_lh_eeg,2)
    slow_gamma_lh_hpc = bandpass_colgin(hpc_lh_eeg(:,z),25,50,250);
    slow_gamma_lh_hpc = hilbert(slow_gamma_lh_hpc);
    hpc_phase_lh(z,:) = angle(slow_gamma_lh_hpc);
    slow_gamma_lh_pfc = bandpass_colgin(pfc_lh_eeg(:,z),5,12,250);
    slow_gamma_lh_pfc = hilbert(slow_gamma_lh_pfc);
    pfc_phase_lh(z,:) = angle(slow_gamma_lh_pfc);
    lh_slow_gamma_ang_diff(z,:) = angdiff(hpc_phase_lh(z,:),pfc_phase_lh(z,:));
end

results.low_high.slow_gamma_phase_hpc = hpc_phase_lh;
results.low_high.slow_gamma_phase_pfc = pfc_phase_lh;
results.low_high.slow_gamma_angle_diff = lh_slow_gamma_ang_diff;
clear hpc_phase_lh pfc_phase_lh lh_slow_gamma_ang_diff
%%%% fast gamma phase
for z = 1:size(hpc_lh_eeg,2)
    fast_gamma_lh_hpc = bandpass_colgin(hpc_lh_eeg(:,z),65,125,250);
    fast_gamma_lh_hpc = hilbert(fast_gamma_lh_hpc);
    hpc_phase_lh(z,:) = angle(fast_gamma_lh_hpc);
    fast_gamma_lh_pfc = bandpass_colgin(pfc_lh_eeg(:,z),5,12,250);
    fast_gamma_lh_pfc = hilbert(fast_gamma_lh_pfc);
    pfc_phase_lh(z,:) = angle(fast_gamma_lh_pfc);
    lh_fast_gamma_ang_diff(z,:) = angdiff(hpc_phase_lh(z,:),pfc_phase_lh(z,:));
end

results.low_high.fast_gamma_phase_hpc = hpc_phase_lh;
results.low_high.fast_gamma_phase_pfc = pfc_phase_lh;
results.low_high.fast_gamma_angle_diff = lh_fast_gamma_ang_diff;
clear hpc_phase_lh pfc_phase_lh lh_fast_gamma_ang_diff

%%% Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hpc_control_eeg = squeeze(hl_eeg.data(1,:,:));
pfc_control_eeg = squeeze(hl_eeg.data(2,:,:));

%%%% theta phase
for z = 1:size(hpc_control_eeg,2)
    theta_control_hpc = bandpass_colgin(hpc_control_eeg(:,z),5,12,250);
    theta_control_hpc = hilbert(theta_control_hpc);
    hpc_phase_control(z,:) = angle(theta_control_hpc);
    theta_control_pfc = bandpass_colgin(pfc_control_eeg(:,z),5,12,250);
    theta_control_pfc = hilbert(theta_control_pfc);
    pfc_phase_control(z,:) = angle(theta_control_pfc);
    control_theta_ang_diff(z,:) = angdiff(hpc_phase_control(z,:),pfc_phase_control(z,:));
end

results.control.theta_phase_hpc = hpc_phase_control;
results.control.theta_phase_pfc = pfc_phase_control;
results.control.theta_angle_diff = control_theta_ang_diff;
clear hpc_phase_control pfc_phase_control control_theta_ang_diff

%%%% slow gamma phase
for z = 1:size(hpc_control_eeg,2)
    slow_gamma_control_hpc = bandpass_colgin(hpc_control_eeg(:,z),25,50,250);
    slow_gamma_control_hpc = hilbert(slow_gamma_control_hpc);
    hpc_phase_control(z,:) = angle(slow_gamma_control_hpc);
    slow_gamma_control_pfc = bandpass_colgin(pfc_control_eeg(:,z),5,12,250);
    slow_gamma_control_pfc = hilbert(slow_gamma_control_pfc);
    pfc_phase_control(z,:) = angle(slow_gamma_control_pfc);
    control_slow_gamma_ang_diff(z,:) = angdiff(hpc_phase_control(z,:),pfc_phase_control(z,:));
end
results.control.slow_gamma_phase_hpc = hpc_phase_control;
results.control.slow_gamma_phase_pfc = pfc_phase_control;
results.control.slow_gamma_angle_diff = control_slow_gamma_ang_diff;
clear hpc_phase_control pfc_phase_control control_slow_gamma_ang_diff

%%%% fast gamma phase
for z = 1:size(hpc_control_eeg,2)
    fast_gamma_control_hpc = bandpass_colgin(hpc_control_eeg(:,z),65,125,250);
    fast_gamma_control_hpc = hilbert(fast_gamma_control_hpc);
    hpc_phase_control(z,:) = angle(fast_gamma_control_hpc);
    fast_gamma_control_pfc = bandpass_colgin(pfc_control_eeg(:,z),5,12,250);
    fast_gamma_control_pfc = hilbert(fast_gamma_control_pfc);
    pfc_phase_control(z,:) = angle(fast_gamma_control_pfc);
    control_fast_gamma_ang_diff(z,:) = angdiff(hpc_phase_control(z,:),pfc_phase_control(z,:));
end

results.control.fast_gamma_phase_hpc = hpc_phase_control;
results.control.fast_gamma_phase_pfc = pfc_phase_control;
results.control.fast_gamma_angle_diff = control_fast_gamma_ang_diff;
clear hpc_phase_control pfc_phase_control control_fast_gamma_ang_diff



%% PHASE




end

