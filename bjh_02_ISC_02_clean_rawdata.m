function [ALLEEG_ISC_clean] = bjh_02_ISC_02_clean_rawdata(ALLEEG_ISC_select_and_reref,params)
%% cleans data using clean_rawdata and interpolates removed channels
% input: ALLEEG_ISC_select_and_reref:    [struct array] structure array containing imported
%                       data from all participants (only first block)
%        params:        [struct] all parameters needed given as subfields
%                       (params.flatline, params.highpass, params.chancrit
%                       params.linenoise, params.burstcrit, params.wincrit,
%                       for detailed description see clean_rawdata.m)
%                       params.visualize: if 1 then vis_artifact is called
% 
% output: ALLEEG_ISC_clean: [struct array] structure array containing cleaned 
%                       data of all participants
% 
% author: Björn Holtze
% date: 14.01.22

for s = 1:size(ALLEEG_ISC_select_and_reref,2)
    
    % high-pass filter (only for comparison before/after with vis_artifacts)
    EEG_filt = clean_drifts(ALLEEG_ISC_select_and_reref(s),params.highpass);
    
    % clean data using default clean_rawdata (default parameters)
    EEG_clean = clean_rawdata(ALLEEG_ISC_select_and_reref(s),params.flatline,params.highpass,params.chancrit,...
        params.linenoise,params.burstcrit,params.wincrit);
    
    if params.visualize
        % compare data before and after cleaning with vis_artifacts
        vis_artifacts(EEG_clean,EEG_filt);
        
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
    
    ALLEEG_ISC_clean(s) = EEG_clean;
    clear EEG_filt;
    clear EEG_clean;
    close;
end

end

