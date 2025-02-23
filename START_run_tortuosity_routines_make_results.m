clear all

if ismac
    cd('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/');
    load('CTL_AND_MIA_DATA.mat');
elseif ispc
   [file] = uigetfile();
   load(file);
end
%cd('/Volumes/BigShare/MIA_Coherence_Bilkey')

for z = 1:size(CTL,1)
    step_one = strsplit(CTL.Filename{z},'_');
    CTL.unique{z} = strcat(step_one{1},step_one{2},step_one{3},step_one{4},step_one{5});
end
[~,ctlvals] = unique(CTL.unique);
CTL_unique = CTL(ctlvals,:);
CTL_unique(13,:) = [];
 
for z = 1:size(MIA,1)
    step_one = strsplit(MIA.Filename{z},'_');
    MIA.unique{z} = strcat(step_one{1},step_one{2},step_one{3},step_one{4},step_one{5});
end
[~,miavals] = unique(MIA.unique);
MIA_unique = MIA(miavals,:);
MIA_unique(4,:) = [];
clear ctlvals miavals step_one
CTL_unique(14,:) = [];

%%% CHOOSE WHICH ANALYSES %%%
lfp_analysis = true;
cohen_analysis = true;
cross_frequency_analysis = true;
phase_analysis = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if lfp_analysis == true
for zz = 1:size(CTL_unique,1)
    ts = CTL_unique.cellTS{zz};
    posx = CTL_unique.posx{zz};
    posy = CTL_unique.posy{zz};
    post = CTL_unique.post{zz};
    eeg_one = CTL_unique.hires_eeg{zz};
    eeg_two = CTL_unique.hires_eeg_2{zz};
    fs = CTL_unique.hires_samplerate(zz);
    Animal = CTL_unique.Animal{zz};
    Date = strsplit(CTL_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(CTL_unique.cellnum(zz)),',Tetrode_',num2str(CTL_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = CTL_unique.lores_eeg{zz};
    lores_eeg_two = CTL_unique.lores_eeg_2{zz};
    lores_fs = CTL_unique.lores_samplerate(zz);
    lores_fs = CTL_unique.lores_samplerate(zz);
    fs = CTL_unique.hires_samplerate(zz);
    [tortuosity_lfp_results_ctl] = do_tortuosity_conjoint_v2(posx,posy,post,ts,lores_eeg_one,lores_eeg_two,animal_id,lores_fs);
    tortuosity_low_high_ctl(zz) = tortuosity_lfp_results_ctl.low_high;
    tortuosity_high_low_ctl(zz) = tortuosity_lfp_results_ctl.high_low;
    tortuosity_control_ctl(zz) = tortuosity_lfp_results_ctl.control;
    tortuosity_freqs.power_freqs = tortuosity_lfp_results_ctl.power_frequencies;
    tortuosity_freqs.transition_timebase = tortuosity_lfp_results_ctl.transition_timebase;
    tortuosity_freqs.coherence_freqs = tortuosity_low_high_ctl.freq_wavelet_coherence_transition;
    fprintf('\ncontrol lfp iteration = %d of %d \n',zz,size(CTL_unique,1));
    clear tortuosity_lfp_results_ctl
end
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/CTL_Tortuosity_LFP_Results_low_high_v3.mat','tortuosity_low_high_ctl');
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/CTL_Tortuosity_LFP_Results_high_low_v3.mat','tortuosity_high_low_ctl');
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/CTL_Tortuosity_LFP_Results_control_v3.mat','tortuosity_control_ctl');
clear tortuosity_low_high_ctl tortuosity_high_low_ctl tortuosity_control_ctl

for zz = 1:size(MIA_unique,1)
    ts = MIA_unique.cellTS{zz};
    posx = MIA_unique.posx{zz};
    posy = MIA_unique.posy{zz};
    post = MIA_unique.post{zz};
    eeg_one = MIA_unique.hires_eeg{zz};
    eeg_two = MIA_unique.hires_eeg_2{zz};
    fs = MIA_unique.hires_samplerate(zz);
    Animal = MIA_unique.Animal{zz};
    Date = strsplit(MIA_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(MIA_unique.cellnum(zz)),',Tetrode_',num2str(MIA_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = MIA_unique.lores_eeg{zz};
    lores_eeg_two = MIA_unique.lores_eeg_2{zz};
    lores_fs = MIA_unique.lores_samplerate(zz);
    lores_fs = MIA_unique.lores_samplerate(zz);
    fs = MIA_unique.hires_samplerate(zz);
    [tortuosity_lfp_results_MIA] = do_tortuosity_conjoint_v2(posx,posy,post,ts,lores_eeg_one,lores_eeg_two,animal_id,lores_fs);
    tortuosity_low_high_MIA(zz) = tortuosity_lfp_results_MIA.low_high;
    tortuosity_high_low_MIA(zz) = tortuosity_lfp_results_MIA.high_low;
    tortuosity_control_MIA(zz) = tortuosity_lfp_results_MIA.control;
    tortuosity_freqs.freqs = tortuosity_lfp_results_MIA.power_frequencies;
    fprintf('\n MIA lfp iteration = %d of %d \n',zz,size(MIA_unique,1));
    clear tortuosity_lfp_results_MIA
end
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/MIA_Tortuosity_LFP_Results_low_high_v3.mat','tortuosity_low_high_MIA');
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/MIA_Tortuosity_LFP_Results_high_low_v3.mat','tortuosity_high_low_MIA');
save('//Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/MIA_Tortuosity_LFP_Results_control_v3.mat','tortuosity_control_MIA');
clear tortuosity_low_high_MIA tortuosity_high_low_MIA tortuosity_control_MIA

for zz = 1:size(CTL_unique,1)
    ts = CTL_unique.cellTS{zz};
    posx = CTL_unique.posx{zz};
    posy = CTL_unique.posy{zz};
    post = CTL_unique.post{zz};
    eeg_one = CTL_unique.hires_eeg{zz};
    eeg_two = CTL_unique.hires_eeg_2{zz};
    fs = CTL_unique.hires_samplerate(zz);
    Animal = CTL_unique.Animal{zz};
    Date = strsplit(CTL_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(CTL_unique.cellnum(zz)),',Tetrode_',num2str(CTL_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = CTL_unique.lores_eeg{zz};
    lores_eeg_two = CTL_unique.lores_eeg_2{zz};
    lores_fs = CTL_unique.lores_samplerate(zz);
    lores_fs = CTL_unique.lores_samplerate(zz);
    fs = CTL_unique.hires_samplerate(zz);
    [spectral_granger_results_ctl(zz)] = do_tortuosity_spectral_granger(posx,posy,post,ts,lores_eeg_one,lores_eeg_two,animal_id,lores_fs);
    fprintf('\ncontrol lfp iteration = %d of %d \n',zz,size(CTL_unique,1));
end
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/CTL_Tortuosity_spectral_granger_Results_v3.mat','spectral_granger_results_ctl')
clear spectral_granger_results_ctl

for zz = 1:size(MIA_unique,1)
    ts = MIA_unique.cellTS{zz};
    posx = MIA_unique.posx{zz};
    posy = MIA_unique.posy{zz};
    post = MIA_unique.post{zz};
    eeg_one = MIA_unique.hires_eeg{zz};
    eeg_two = MIA_unique.hires_eeg_2{zz};
    fs = MIA_unique.hires_samplerate(zz);
    Animal = MIA_unique.Animal{zz};
    Date = strsplit(MIA_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(MIA_unique.cellnum(zz)),',Tetrode_',num2str(MIA_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = MIA_unique.lores_eeg{zz};
    lores_eeg_two = MIA_unique.lores_eeg_2{zz};
    lores_fs = MIA_unique.lores_samplerate(zz);
    lores_fs = MIA_unique.lores_samplerate(zz);
    fs = MIA_unique.hires_samplerate(zz);
    [spectral_granger_results_mia(zz)] = do_tortuosity_spectral_granger(posx,posy,post,ts,lores_eeg_one,lores_eeg_two,animal_id,lores_fs);
    fprintf('\nmia lfp iteration = %d of %d \n',zz,size(MIA_unique,1));
end
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/MIA_Tortuosity_spectral_granger_Results_v3.mat','spectral_granger_results_mia')
clear spectral_granger_results_mia
else 
end

if cohen_analyses == true
% Cohen analyses
for zz = 1:size(CTL_unique,1)
    ts = CTL_unique.cellTS{zz};
    posx = CTL_unique.posx{zz};
    posy = CTL_unique.posy{zz};
    post = CTL_unique.post{zz};
    eeg_one = CTL_unique.hires_eeg{zz};
    eeg_two = CTL_unique.hires_eeg_2{zz};
    fs = CTL_unique.hires_samplerate(zz);
    Animal = CTL_unique.Animal{zz};
    Date = strsplit(CTL_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(CTL_unique.cellnum(zz)),',Tetrode_',num2str(CTL_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = CTL_unique.lores_eeg{zz};
    lores_eeg_two = CTL_unique.lores_eeg_2{zz};
    lores_fs = CTL_unique.lores_samplerate(zz);
    lores_fs = CTL_unique.lores_samplerate(zz);
    fs = CTL_unique.hires_samplerate(zz);
    plotfig = 0;
    %%%%% LOW TO HIGH TRANSITIONS
    [EEG_low_high] = get_transitions_low_high(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [erp_results_low_high_CTL(zz)] = calculate_ERP_phase_locking(EEG_low_high,animal_id,plotfig);
    [cohen_cohere_results_low_high_CTL(zz)] = calculate_cohen_spec_coherence(EEG_low_high,animal_id,plotfig);
    [cohen_pli_results_low_high_CTL(zz)] = calculate_phase_lag_index(EEG_low_high,animal_id,plotfig);
    [cohen_power_correlation_results_low_high_CTL(zz)] = calculate_power_correlation(EEG_low_high,animal_id,plotfig);

    %%%%% HIGH TO LOW TRANSITIONS
    [EEG_high_low] = get_transitions_high_low(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [erp_results_high_low_CTL(zz)] = calculate_ERP_phase_locking(EEG_high_low,animal_id,plotfig);
    [cohen_cohere_results_high_low_CTL(zz)] = calculate_cohen_spec_coherence(EEG_high_low,animal_id,plotfig);
    [cohen_pli_results_high_low_CTL(zz)] = calculate_phase_lag_index(EEG_high_low,animal_id,plotfig);
    [cohen_power_correlation_results_high_low_CTL(zz)] = calculate_power_correlation(EEG_high_low,animal_id,plotfig);

    fprintf('\n CTL lfp iteration = %d of %d \n',zz,size(CTL_unique,1));
end

save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/CTL_master_cohen_results_v3.mat','erp_results_low_high_CTL','erp_results_high_low_CTL'...
    ,'cohen_cohere_results_low_high_CTL','cohen_cohere_results_high_low_CTL','cohen_pli_results_high_low_CTL','cohen_pli_results_low_high_CTL'...
    ,'cohen_power_correlation_results_low_high_CTL','cohen_power_correlation_results_high_low_CTL');

%%% MIA COHEN
% Event related potential results, Phase lag results, coherence results, power correlation
% results
for zz = 1:size(MIA_unique,1)
    ts = MIA_unique.cellTS{zz};
    posx = MIA_unique.posx{zz};
    posy = MIA_unique.posy{zz};
    post = MIA_unique.post{zz};
    eeg_one = MIA_unique.hires_eeg{zz};
    eeg_two = MIA_unique.hires_eeg_2{zz};
    fs = MIA_unique.hires_samplerate(zz);
    Animal = MIA_unique.Animal{zz};
    Date = strsplit(MIA_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(MIA_unique.cellnum(zz)),',Tetrode_',num2str(MIA_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = MIA_unique.lores_eeg{zz};
    lores_eeg_two = MIA_unique.lores_eeg_2{zz};
    lores_fs = MIA_unique.lores_samplerate(zz);
    lores_fs = MIA_unique.lores_samplerate(zz);
    fs = MIA_unique.hires_samplerate(zz);
    plotfig = 0;
    %%%%% LOW TO HIGH TRANSITIONS
    [EEG_low_high] = get_transitions_low_high(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [erp_results_low_high_mia(zz)] = calculate_ERP_phase_locking(EEG_low_high,animal_id,plotfig);
    [cohen_cohere_results_low_high_mia(zz)] = calculate_cohen_spec_coherence(EEG_low_high,animal_id,plotfig);
    [cohen_pli_results_low_high_mia(zz)] = calculate_phase_lag_index(EEG_low_high,animal_id,plotfig);
    [cohen_power_correlation_results_low_high_mia(zz)] = calculate_power_correlation(EEG_low_high,animal_id,plotfig);
    %%%%% HIGH TO LOW TRANSITIONS
    [EEG_high_low] = get_transitions_high_low(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [erp_results_high_low_mia(zz)] = calculate_ERP_phase_locking(EEG_high_low,animal_id,plotfig);
    [cohen_cohere_results_high_low_mia(zz)] = calculate_cohen_spec_coherence(EEG_high_low,animal_id,plotfig);
    [cohen_pli_results_high_low_mia(zz)] = calculate_phase_lag_index(EEG_high_low,animal_id,plotfig);
    [cohen_power_correlation_results_high_low_mia(zz)] = calculate_power_correlation(EEG_high_low,animal_id,plotfig);
    fprintf('\n MIA lfp iteration = %d of %d \n',zz,size(MIA_unique,1));
end

save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/MIA_master_cohen_results_v3.mat','erp_results_low_high_mia','erp_results_high_low_mia'...
    ,'cohen_cohere_results_low_high_mia','cohen_cohere_results_high_low_mia','cohen_pli_results_high_low_mia','cohen_pli_results_low_high_mia'...
    ,'cohen_power_correlation_results_low_high_mia','cohen_power_correlation_results_high_low_mia');

else
end


if cross_freq_analysis == true

%%%% CROSS FREQUENCY COUPLING
for zz = 1:size(CTL_unique,1)
    ts = CTL_unique.cellTS{zz};
    posx = CTL_unique.posx{zz};
    posy = CTL_unique.posy{zz};
    post = CTL_unique.post{zz};
    eeg_one = CTL_unique.hires_eeg{zz};
    eeg_two = CTL_unique.hires_eeg_2{zz};
    fs = CTL_unique.hires_samplerate(zz);
    Animal = CTL_unique.Animal{zz};
    Date = strsplit(CTL_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(CTL_unique.cellnum(zz)),',Tetrode_',num2str(CTL_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = CTL_unique.lores_eeg{zz};
    lores_eeg_two = CTL_unique.lores_eeg_2{zz};
    lores_fs = CTL_unique.lores_samplerate(zz);
    lores_fs = CTL_unique.lores_samplerate(zz);
    fs = CTL_unique.hires_samplerate(zz);
    plotfig = 0;
    [EEG_low_high] = get_transitions_low_high(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [cross_freq_results_low_high_CTL(zz)] = calculate_cross_frequency_coupling(EEG_low_high,animal_id,plotfig);
    [EEG_high_low] = get_transitions_high_low(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [cross_freq_results_high_low_CTL(zz)] = calculate_cross_frequency_coupling(EEG_high_low,animal_id,plotfig);
        fprintf('\n CTL lfp iteration = %d of %d \n',zz,size(CTL_unique,1));
end
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/CTL_master_cross_freq_results_v3.mat','cross_freq_results_low_high_CTL','cross_freq_results_high_low_CTL')


% MIA


for zz = 1:size(MIA_unique,1)
    ts = MIA_unique.cellTS{zz};
    posx = MIA_unique.posx{zz};
    posy = MIA_unique.posy{zz};
    post = MIA_unique.post{zz};
    eeg_one = MIA_unique.hires_eeg{zz};
    eeg_two = MIA_unique.hires_eeg_2{zz};
    fs = MIA_unique.hires_samplerate(zz);
    Animal = MIA_unique.Animal{zz};
    Date = strsplit(MIA_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(MIA_unique.cellnum(zz)),',Tetrode_',num2str(MIA_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = MIA_unique.lores_eeg{zz};
    lores_eeg_two = MIA_unique.lores_eeg_2{zz};
    lores_fs = MIA_unique.lores_samplerate(zz);
    lores_fs = MIA_unique.lores_samplerate(zz);
    fs = MIA_unique.hires_samplerate(zz);
    plotfig = 0;
    [EEG_low_high] = get_transitions_low_high(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [cross_freq_results_low_high_MIA(zz)] = calculate_cross_frequency_coupling(EEG_low_high,animal_id,plotfig);
    [EEG_high_low] = get_transitions_high_low(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [cross_freq_results_high_low_MIA(zz)] = calculate_cross_frequency_coupling(EEG_high_low,animal_id,plotfig);
        fprintf('\n MIA lfp iteration = %d of %d \n',zz,size(MIA_unique,1));
end
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/MIA_master_cross_freq_results_v3.mat','cross_freq_results_low_high_MIA','cross_freq_results_high_low_MIA')

else
end


%%%%%%% PHASE ANALYSIS %%%%%%%%%%%%

% CTL
if phase_analysis == true

for zz = 1:size(CTL_unique,1)
    ts = CTL_unique.cellTS{zz};
    posx = CTL_unique.posx{zz};
    posy = CTL_unique.posy{zz};
    post = CTL_unique.post{zz};
    eeg_one = CTL_unique.hires_eeg{zz};
    eeg_two = CTL_unique.hires_eeg_2{zz};
    fs = CTL_unique.hires_samplerate(zz);
    Animal = CTL_unique.Animal{zz};
    Date = strsplit(CTL_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(CTL_unique.cellnum(zz)),',Tetrode_',num2str(CTL_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = CTL_unique.lores_eeg{zz};
    lores_eeg_two = CTL_unique.lores_eeg_2{zz};
    lores_fs = CTL_unique.lores_samplerate(zz);
    lores_fs = CTL_unique.lores_samplerate(zz);
    fs = CTL_unique.hires_samplerate(zz);
    plotfig = 0;
    [EEG_low_high] = get_transitions_low_high(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [EEG_high_low] = get_transitions_high_low(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [phase_results_CTL(zz)] = calculate_phase_and_phase_slope_index(EEG_low_high,EEG_high_low);
        fprintf('\n CTL lfp iteration = %d of %d \n',zz,size(CTL_unique,1));
end    
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/CTL_master_phase_results_v3.mat','phase_results_CTL')
%MIA
for zz = 1:size(MIA_unique,1)
    ts = MIA_unique.cellTS{zz};
    posx = MIA_unique.posx{zz};
    posy = MIA_unique.posy{zz};
    post = MIA_unique.post{zz};
    eeg_one = MIA_unique.hires_eeg{zz};
    eeg_two = MIA_unique.hires_eeg_2{zz};
    fs = MIA_unique.hires_samplerate(zz);
    Animal = MIA_unique.Animal{zz};
    Date = strsplit(MIA_unique.Filename{zz},'SP');
    Date = Date{1};
    cell_tet = strcat('Cell_',num2str(MIA_unique.cellnum(zz)),',Tetrode_',num2str(MIA_unique.tetrode(zz)));
    animal_id = strcat(Animal,'_',Date,cell_tet);
    lores_eeg_one = MIA_unique.lores_eeg{zz};
    lores_eeg_two = MIA_unique.lores_eeg_2{zz};
    lores_fs = MIA_unique.lores_samplerate(zz);
    lores_fs = MIA_unique.lores_samplerate(zz);
    fs = MIA_unique.hires_samplerate(zz);
    plotfig = 0;
    [EEG_low_high] = get_transitions_low_high(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [EEG_high_low] = get_transitions_high_low(posx,posy,lores_eeg_one,lores_eeg_two,lores_fs);
    [phase_results_MIA(zz)] = calculate_phase_and_phase_slope_index(EEG_low_high,EEG_high_low);
        fprintf('\n MIA lfp iteration = %d of %d \n',zz,size(MIA_unique,1));
end    
save('/Volumes/Drivepool/Rob/Science/Active_Projects/Tortuosity/Results/MIA_master_phase_results_v3.mat','phase_results_MIA')