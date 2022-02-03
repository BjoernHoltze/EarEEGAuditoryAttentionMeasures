function bjh_01_AAD_04_extract_envelope_mirkovic2016(file_in_audio,DATAPATH,params)
%% extract speech envelope from the .wav files of each block
% First, the absolute Hilbert transform of the .wav file is computed, then
% it is lowpass filtered, then it is normalized and downsampled
% 
% input:    files_in_audiobook:     struct created with dir() which lists all .wav files
%                                   to be loaded. From these .wav files the
%                                   speech envelope is extracted.
%           path_out:               [string] path where .mat files containing the 
%                                   extracted envelopes will be stored
%           params:                 structure containing field lp_filter
%                                   and new sampling rate
% 
% author: Björn Holtze
% date: 16.04.2021

for f = 1:size(file_in_audio,1)
    disp(['Extracting envelope from ',file_in_audio(f).name,' ...']);
    
    % load raw audio signal
    [wav, fs_old] = audioread([file_in_audio(f).folder,filesep,...
        file_in_audio(f).name]);    
    
    % normalize
    wav_norm = wav/std(wav);
    
    % compute absolute Hilbert transform
    abs_hilb = abs(hilbert(wav_norm));
    
    % lowpass filter absolute Hilbert transform
    f_order = 3;  
    Wn = params.lp_filter/(fs_old/2);
    [b1,a1] = butter(f_order,Wn,'low');
    env = filtfilt(b1,a1,abs_hilb);
    
    % downsample
    env_resampled = resample(env,params.fs_new,fs_old);
    
    % visualize each step (only consider first minute)
    if params.visualize_envelope
        vis_win = 3; % in s
        figure('Name',file_in_audio(f).name,'Units','centimeters',...
            'Position',[0.5,2,12,16]);
        subplot(5,1,1);
        plot(linspace(0,1,vis_win*fs_old),wav(1:vis_win*fs_old));
        title('Original Speech');
        subplot(5,1,2);
        plot(linspace(0,1,vis_win*fs_old),wav_norm(1:vis_win*fs_old));
        title('Normalized');
        subplot(5,1,3);
        plot(linspace(0,1,vis_win*fs_old),abs_hilb(1:vis_win*fs_old));
        title('Absolute Hilbert Transform');
        subplot(5,1,4);
        plot(linspace(0,1,vis_win*fs_old),env(1:vis_win*fs_old));
        title(['Lowpass Filtered (Butterworth: ',num2str(params.lp_filter),...
            ' Hz, order: ',num2str(f_order),')']);
        subplot(5,1,5);
        plot(linspace(0,1,vis_win*params.fs_new),env_resampled(1:vis_win*params.fs_new));
        title('Resampled');
    end
    
    % pause for a second so that figure is plotted
    pause(1);
    
    % save processed envelope
    save([DATAPATH,'Envelope_Mirkovic2016',file_in_audio(f).name(7:10),'.mat'],'env_resampled');
end
    
    
end

