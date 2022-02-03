function [ALLEEG_SpE_clean] = bjh_03_SpE_02_clean_rawdata(ALLEEG_SpE_select_and_reref,params)
%% cleans data using clean_rawdata and interpolates removed channels
% 
% input: ALLEEG_SpE_select_and_reref:    
%                       [struct array] structure array containing imported
%                       data from all participants (all relevant blocks
%                       concatenated)
%           params:     [struct] all parameters needed given as subfields
%                       (params.flatline, params.highpass, params.chancrit
%                       params.linenoise, params.burstcrit, params.wincrit,
%                       for detailed description see clean_rawdata.m)
%                       params.visualize: if 1 then vis_artifact is called
% 
% output: ALLEEG_SpE_clean: [struct array] structure array containing cleaned 
%                       data of all participants
% 
% author: Björn Holtze
% date: 14.01.22

for s = 1:size(ALLEEG_SpE_select_and_reref,2)
    
    % high-pass filter (only for comparison before/after with vis_artifacts)
    EEG_filt = clean_drifts(ALLEEG_SpE_select_and_reref(s),params.highpass);
    
    % identify channels only weakly correlate with other channels
    [~, bad_chans] = clean_channels(EEG_filt,params.chancrit,100);
    
    % clean data using default clean_rawdata (default parameters)
    EEG_clean = clean_rawdata(ALLEEG_SpE_select_and_reref(s),params.flatline,params.highpass,'off',...
        params.linenoise,params.burstcrit,params.wincrit);
    
    % write removed channels into EEG_clean
    EEG_clean.bad_chans = {EEG_clean.chanlocs(bad_chans).labels};
    
    if params.visualize
        % compare data before and after cleaning with vis_artifacts
        vis_artifacts(EEG_filt,EEG_clean);
        
        % vis_artifacts needs some time to load figure
        pause(1);
        
        while 1
            i = input('Enter "1" to continue: ');
            if i == 1
                break;
            end
        end
    end
    
    EEG_clean.params = params;
    
    ALLEEG_SpE_clean(s) = EEG_clean;
    clear EEG_filt;
    clear EEG_clean;
    close;
end
end

