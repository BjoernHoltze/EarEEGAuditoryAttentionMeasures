function [ALLEEG_ISC_filt] = bjh_02_ISC_03_filt_downsample(ALLEEG_ISC_clean,params)
%% filters, normalizes and downsamples the cEEGrid EEG 
% input:    ALLEEG_ISC_clean:  structure array containing EEG structures from each
%                       partcipant (similar to ALLEEG in eeglab), cleaned with ASR
%           params:     structure containing highpass and lowpass filter
%                       edges as well as new sampling rate as subfields
% 
% output:   [ALLEEG_ISC_filt: structure array containing filtered, normalized, 
%                        and downsampled EEG structures from all subjects
%
% author: Björn Holtze
% date: 14.01.2022

for s = 1:size(ALLEEG_ISC_clean,2)
    disp(['Processing subject ',num2str(s),' ...']);
    
    % lowpass filter data
    EEG = pop_firws(ALLEEG_ISC_clean(s),'fcutoff',params.lp_filter,'ftype','lowpass',...
        'wtype','hann','forder',100);
    
    % resample
    EEG = pop_resample(EEG,params.fs_new);
    
    % highpass filter data
    EEG = pop_firws(EEG,'fcutoff',params.hp_filter,'ftype','highpass',...
        'wtype','hann','forder',500);
    
    % epoch according to blocks (0 to 600 s)
    EEG = pop_epoch(EEG,{'StartTrigger'},[0  600],'epochinfo','yes');
    
    % concatenate epochs
    EEG = eeg_epoch2continuous(EEG);
    
    % note parameters
    EEG.lp_filter = params.lp_filter;
    EEG.hp_filter = params.hp_filter;
    
    ALLEEG_ISC_filt(s) = EEG;
end

end

