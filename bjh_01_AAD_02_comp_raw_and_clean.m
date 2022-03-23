function [PSD_change_db] = bjh_01_AAD_02_comp_raw_and_clean(fileout_AAD_psd_raw_vs_clean,...
    ALLEEG_AAD_select_and_reref,ALLEEG_AAD_clean,seg_s,freq_range)

%% computes PSD change from before to after artifact correction
% input:    fileout_AAD_spec_pow_change_plot_supp_figure:      
%               full path including filename of figure to be stored
%           ALLEEG_AAD_select_and_reref:    ALLEEG structure before artifact correction
%           ALLEEG_AAD_clean:               ALLEEG structure after artifact correction
%           seg_s:                          length of segments (in seconds)
%           freq_range:                     frequency band to be considered (2-8 Hz)
% 
% output:   PSD_change_db:                  (subjects x segments) 
%                                           PSD change in dB in the specified freuency band
%                                           from before to after artifact correction
%                                           averaged over channels (with participants)
% 
% author: Björn Holtze 
% date: 18.03.2022

for s = 1:size(ALLEEG_AAD_select_and_reref,2)
    % remove periods before start and after end of 10-minute block
    EEG_raw = pop_epoch(ALLEEG_AAD_select_and_reref(s),{'StartTrigger'},[0  600],'epochinfo','yes');
    EEG_clean = pop_epoch(ALLEEG_AAD_clean(s),{'StartTrigger'},[0  600],'epochinfo','yes');
    
    % concatenate epochs
    EEG_raw = eeg_epoch2continuous(EEG_raw);
    EEG_clean = eeg_epoch2continuous(EEG_clean);
    
    % remove boundary events
    EEG_raw.event = EEG_raw.event(~strcmp({EEG_raw.event.type},'boundary'));
    EEG_clean.event = EEG_clean.event(~strcmp({EEG_clean.event.type},'boundary'));
    
    % filter data
    EEG_raw = pop_firws(EEG_raw,'fcutoff',freq_range(2),'ftype','lowpass',...
        'wtype','hann','forder',100);
    EEG_raw = pop_firws(EEG_raw,'fcutoff',freq_range(1),'ftype','highpass',...
        'wtype','hann','forder',500);
    EEG_clean = pop_firws(EEG_clean,'fcutoff',freq_range(2),'ftype','lowpass',...
        'wtype','hann','forder',100);
    EEG_clean = pop_firws(EEG_clean,'fcutoff',freq_range(1),'ftype','highpass',...
        'wtype','hann','forder',500);
    
    % split data into continuous segments
    ALLEEG_raw_ep(s) = eeg_regepochs(EEG_raw,'recurrence',seg_s,'limits',[0,seg_s]);
    ALLEEG_clean_ep(s) = eeg_regepochs(EEG_clean,'recurrence',seg_s,'limits',[0,seg_s]);
    
    % compute PSD with FFT (raw)
    N_raw = size(ALLEEG_raw_ep(s).data,2);
    Fs_raw = ALLEEG_raw_ep(s).srate;
    xdft_raw_mirror(s,:,:,:) = fft(ALLEEG_raw_ep(s).data,N_raw,2);
    xdft_raw(s,:,:,:) = xdft_raw_mirror(s,:,1:N_raw/2+1,:);
    psdx_raw = (1/(Fs_raw*N_raw)) * abs(xdft_raw).^2;
    psdx_raw(s,:,2:end-1,:) = 2*psdx_raw(s,:,2:end-1,:);
    
    % compute PSD with FFT (clean)
    N_clean = size(ALLEEG_clean_ep(s).data,2);
    Fs_clean = ALLEEG_clean_ep(s).srate;
    xdft_clean_mirror(s,:,:,:) = fft(ALLEEG_clean_ep(s).data,N_clean,2);
    xdft_clean(s,:,:,:) = xdft_clean_mirror(s,:,1:N_clean/2+1,:);
    psdx_clean = (1/(Fs_clean*N_clean)) * abs(xdft_clean).^2;
    psdx_clean(s,:,2:end-1,:) = 2*psdx_clean(s,:,2:end-1,:);
    
    % compute and plot histogram with power change in dB
    freq_vec = 0:Fs_raw/N_raw:Fs_raw/2;
    freq_vec_range = freq_vec >= freq_range(1) & freq_vec <= freq_range(2);
    spec_pow_raw(s,:) = sum(squeeze(mean(psdx_raw(s,:,freq_vec_range,:),2)),1); % average over channels
    spec_pow_clean(s,:) = sum(squeeze(mean(psdx_clean(s,:,freq_vec_range,:),2)),1); % average over channels
    PSD_change_db(s,:) = 10*log10(spec_pow_clean(s,:)./spec_pow_raw(s,:));
    
    clear EEG_raw EEG_clean
end

    figure('units','centimeters','outerposition',[5 5 18 10]);
    histogram(PSD_change_db(:),-7:0.1:2,'Normalization','probability');
    xlabel('Spectral Power Change [dB]');
    ylabel('Proportion of 1-Second Segments [%]');
    title(['Spectral Power Change (2-8 Hz)',newline,'From Before to After Artifact Correction']);
    box off;
    ytix = get(gca, 'YTick');
    set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
    
    saveas(gcf,[fileout_AAD_psd_raw_vs_clean,'.jpg']);
    saveas(gcf,[fileout_AAD_psd_raw_vs_clean,'.png']);
    saveas(gcf,[fileout_AAD_psd_raw_vs_clean,'.svg']);
    saveas(gcf,[fileout_AAD_psd_raw_vs_clean,'.eps'],'epsc');
    saveas(gcf,[fileout_AAD_psd_raw_vs_clean,'.tiff']);
    saveas(gcf,[fileout_AAD_psd_raw_vs_clean,'.pdf']);
    close;
    
end

