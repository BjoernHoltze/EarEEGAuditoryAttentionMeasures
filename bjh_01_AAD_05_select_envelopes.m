function [ALLEEG_AAD_env] = bjh_01_AAD_05_select_envelopes(file_in_conc_env,ALLEEG_prep)
%% concatenates subject-specifically selected blocks of the envelope
% 
% input:    file_in_env:    struct created with dir.m listing all processed
%                           envelopes (2 streams 5 blocks each)
%           ALLEEG_prep:    struct containing all preprocessed EEG datasets
%                           (similar to ALLEEG structure in EEGLAB)
% 
% output:   ALLEEG_env:     struct containing all information that ALLEEG_prep 
%                           contained plus the selected blocks of the envelop 
%                           concatenated as subfield to ALLEEG_env 
%                           ALLEEG_env.env_att
%                           ALLEEG_env.env_ign
%
% author: Björn Holtze 
% date: 21.04.2021

% preallocate space (600 seconds x sampling rate, 5 blocks)
stream_1 = zeros(600*ALLEEG_prep(1).srate,5);
stream_2 = zeros(600*ALLEEG_prep(1).srate,5);

% there are 10 files in file_in_env (2 stories with 5 blocks each)
for f = 1:size(file_in_conc_env,1)
    % load extracted envelope 
    load([file_in_conc_env(f).folder,filesep,file_in_conc_env(f).name],'env_resampled');
    
    % only consider exactly 10 min per block
    if strcmp(file_in_conc_env(f).name(end-6),'1')
        stream_1(:,str2double(file_in_conc_env(f).name(end-4))) = ...
            env_resampled(1:600*ALLEEG_prep(1).srate,1);
    elseif strcmp(file_in_conc_env(f).name(end-6),'2')
        stream_2(:,str2double(file_in_conc_env(f).name(end-4))) = ...
            env_resampled(1:600*ALLEEG_prep(1).srate,1);
    end
end

ALLEEG_AAD_env = ALLEEG_prep;

for s = 1:size(ALLEEG_prep,2)
    
    if ALLEEG_prep(s).attended_ch == 1
        % concatenate subjet specific blocks of attended stream
        ALLEEG_AAD_env(s).env_att = reshape(stream_1(:,ALLEEG_AAD_env(s).selected_bl),...
            [1,600*ALLEEG_prep(1).srate*3]);
        % concatenate subjet specific blocks of ignored stream
        ALLEEG_AAD_env(s).env_ign = reshape(stream_2(:,ALLEEG_AAD_env(s).selected_bl),...
            [1,600*ALLEEG_prep(1).srate*3]);
    elseif ALLEEG_prep(s).attended_ch == 2
        % concatenate subjet specific blocks of attended stream
        ALLEEG_AAD_env(s).env_att = reshape(stream_2(:,ALLEEG_AAD_env(s).selected_bl),...
            [1,600*ALLEEG_prep(1).srate*3]);
        % concatenate subjet specific blocks of ignored stream
        ALLEEG_AAD_env(s).env_ign = reshape(stream_1(:,ALLEEG_AAD_env(s).selected_bl),...
            [1,600*ALLEEG_prep(1).srate*3]);
    end
end

end

