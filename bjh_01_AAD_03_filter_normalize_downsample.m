function [ALLEEG_prep] = bjh_01_AAD_03_filter_normalize_downsample(ALLEEG_in,params)
%% filters, normalizes and downsamples the cEEGrid EEG 
% 
% input:    ALLEEG_in:  structure array containing EEG structures from each
%                       partcipant (similar to ALLEEG in eeglab)
%           file_out:   filename (without file extension including full 
%                       path where ALLEEG_prep will be stored
%           params:     structure containing highpass and lowpass filter
%                       edges as well as new sampling rate as subfields
% 
% output:   ALLEEG_prep: structure array containing filetered, normalized, 
%                        and downsampled EEG structures from all subjects
%
% author: Björn Holtze
% date: 16.04.2021

for s = 1:size(ALLEEG_in,2)
    disp(['Processing subject ',num2str(s),' ...']);
    
    % lowpass filter data
    EEG = pop_firws(ALLEEG_in(s),'fcutoff',params.lp_filter,'ftype','lowpass',...
        'wtype','hann','forder',100);
    
    % highpass filter data
    EEG = pop_firws(EEG,'fcutoff',params.hp_filter,'ftype','highpass',...
        'wtype','hann','forder',500);
    
    % Normalize (dividing by standard deviation)
    EEG.data = EEG.data/std(EEG.data(:));
    
    % resample
    EEG = pop_resample(EEG,params.fs_new);
    
    % epoch according to blocks (0 to 600 s)
    EEG = pop_epoch(EEG,{'StartTrigger'},[0  600],'epochinfo','yes');
    
    % concatenate epochs
    EEG = eeg_epoch2continuous(EEG);
    
    % note parameters
    EEG.lp_filter = params.lp_filter;
    EEG.hp_filter = params.hp_filter;
    
    ALLEEG_prep(s) = EEG;
end
end

