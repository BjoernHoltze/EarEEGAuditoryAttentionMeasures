function [SpE_results] = bjh_03_SpE_04_compute_entropy(ALLEEG_SpE_filt,params)
%% computes spectra and their entropy over time for different EEG channels
% input:    ALLEEG_SpE_filt:   structure containing ASR cleaned data
%           params:                 structure containing parameters to
%                                   compute spectral entropy
%               params.entr_freq_min:   minimal frequency band to be considered
%               params.entr_freq_max:   maximal frequency band to be considered
%               params.seg_length:      length of segments in which entropy
%                                       should be computed (in seconds)
% 
% output:   SpE_results: 
%               results.freq_vec:       frequency vector of computed spectrum 
%               results.psd_norm:       normed power spectral density
%                                       (subj x seg x freq x chan)
%               results.entropy_norm:   normed spectral entropy 
%                                       (subj x seg x chan)
%
% author: Björn Holtze (based on Lesenfants et al. 2020)
% date: 14.02.22


fs = ALLEEG_SpE_filt(1).srate;
nb_ch = ALLEEG_SpE_filt(1).nbchan;
nb_subj = size(ALLEEG_SpE_filt,2);
seg_start = 1:fs*params.seg_length:ALLEEG_SpE_filt(1).pnts;
nb_seg = size(seg_start,2);
entropy_norm = zeros(nb_subj,nb_seg,nb_ch);
psd = zeros(nb_subj,nb_seg,fs+1,nb_ch);
psd_norm = zeros(size(psd));
psd_norm_chan_avg = zeros(nb_subj,nb_seg,fs+1);
entropy_norm_chan_avg = zeros(nb_subj,nb_seg);

for subj = 1:nb_subj
    disp(['Processing subject ',num2str(subj),' ...']);
    
    for seg = 1:nb_seg
        % compute multitaper power spectral density (PSD)
        ntaper = 7;
        nw = (ntaper+1)/2;
        nfft = fs;
        [psd(subj,seg,:,:),freq_vec] = pmtm(ALLEEG_SpE_filt(subj).data(:,...
            seg_start(seg):seg_start(seg)+params.seg_length*fs-1)',nw,nfft*2,fs);
        
        entropy = zeros(nb_ch,1);
        for c = 1:nb_ch
            % normalize PSD (divide each frequency by sum over all frequencies)
            psd_norm(subj,seg,:,c) = psd(subj,seg,:,c)/sum(psd(subj,seg,...
                freq_vec >= params.entr_freq_min & ...
                freq_vec <= params.entr_freq_max,c));
            
            % compute entropy
            freq_range =  params.entr_freq_min:0.5:params.entr_freq_max;
            for f = 1:size(freq_range,2)
                entropy(c) = entropy(c) + psd_norm(subj,seg,freq_vec == freq_range(f),c) * ...
                    log(1/psd_norm(subj,seg,freq_vec == freq_range(f),c));
            end
        end
        % normalize entropy by number of frequencies
        entropy_norm(subj,seg,:) = entropy/log(size(freq_range,2))';
        
        % compute PDS averaged over channels 
        psd_norm_chan_avg(subj,seg,:) = mean(psd_norm(subj,seg,:,...
            ~ismember({ALLEEG_SpE_filt(subj).chanlocs.labels},...
            ALLEEG_SpE_filt(subj).bad_chans)),4);
        
        % compute entropy of PSD averaged over channels
        entropy_chan_avg = 0;
        for f = 1:size(freq_range,2)
            entropy_chan_avg = entropy_chan_avg + psd_norm_chan_avg(subj,seg,freq_vec == freq_range(f)) * ...
                log(1/psd_norm_chan_avg(subj,seg,freq_vec == freq_range(f)));
        end    
        
        % normalize entropy over all channels
        entropy_norm_chan_avg(subj,seg) = entropy_chan_avg/log(size(freq_range,2));
    end
end

SpE_results.freq_vec = freq_vec;
SpE_results.psd_norm = psd_norm; 
SpE_results.entropy_norm = entropy_norm;
SpE_results.psd_norm_chan_avg = psd_norm_chan_avg;
SpE_results.entropy_norm_chan_avg = entropy_norm_chan_avg;
for s = 1: nb_subj
    SpE_results.subj_info(s).setname = ALLEEG_SpE_filt(s).setname(1:7);
    SpE_results.subj_info(s).selected_bl = ALLEEG_SpE_filt(s).selected_bl;
    SpE_results.subj_info(s).bad_chans = ALLEEG_SpE_filt(s).bad_chans;
    SpE_results.subj_info(s).chanlocs = ALLEEG_SpE_filt(s).chanlocs;
end
end



