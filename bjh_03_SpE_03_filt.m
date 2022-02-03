function [ALLEEG_SpE_filt] = bjh_03_SpE_03_filt(ALLEEG_SpE_clean,params)
%% filters, normalizes and downsamples the cEEGrid EEG 
% 
% input:    ALLEEG_SpE_clean:  structure array containing EEG structures from each
%                       partcipant (similar to ALLEEG in eeglab)
%           params:     structure containing highpass and lowpass filter
%                       edges as well as new sampling rate as subfields
% 
% output:   ALLEEG_SpE_filt: structure array containing filetered, normalized, 
%                        and downsampled EEG structures from all subjects
%
% author: Björn Holtze
% date: 14.02.2022

for s = 1:size(ALLEEG_SpE_clean,2)
    disp(['Processing subject ',num2str(s),' ...']);
    
    % lowpass filter data
    EEG = pop_firws(ALLEEG_SpE_clean(s),'fcutoff',params.lp_filter,'ftype','lowpass',...
        'wtype','hann','forder',100);
    
    % highpass filter data
    EEG = pop_firws(EEG,'fcutoff',params.hp_filter,'ftype','highpass',...
        'wtype','hann','forder',776);
    
    % epoch according to blocks (0 to 600 s)
    EEG = pop_epoch(EEG,{'StartTrigger'},[0  600],'epochinfo','yes');
    
    % concatenate epochs
    EEG = eeg_epoch2continuous(EEG);
    
    % note parameters
    EEG.lp_filter = params.lp_filter;
    EEG.hp_filter = params.hp_filter;
    
    ALLEEG_SpE_filt(s) = EEG;
end

end

