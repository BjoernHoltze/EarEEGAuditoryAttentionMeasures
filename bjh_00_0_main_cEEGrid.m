%% Analysis Pipeline - Holtze et al. (2022)
% author: Björn Holtze
% date: 14.01.2022

MAINPATH = [pwd,filesep,'..',filesep];
BIDS_folder_name = 'bjh_cEEGrid_attention';
CONFIGPATH = [MAINPATH,BIDS_folder_name,filesep,'sourcedata',filesep];

% add EEGLAB to MATLAB path
addpath(fullfile([MAINPATH,'software',filesep,'eeglab2020_0',filesep]));
% add cEEGrid plugin
addpath(fullfile([MAINPATH,'software',filesep,'eeglab2020_0',filesep,'plugins',filesep,...
    'ceegridplugin-master',filesep,'ceegridplugin-master',filesep]));
% open and close EEGLAB so that functions are included in path as well
eeglab; close 'EEGLAB v2020.0';
% add mTRF Toolbox
addpath(genpath([MAINPATH,'software',filesep,'mTRF-Toolbox-2.1',filesep,...
    'mickcrosse-mTRF-Toolbox-a2806ca',filesep]));

% indicate which steps to skip (which steps have been computed before)
config.AAD_select_and_reref = 1;
config.AAD_clean = 1;
config.AAD_prep = 1;
config.AAD_envelope  = 1;
config.AAD_select_envelope = 1;
config.AAD_val_test_sep = 1;
config.AAD_val_test_mix = 1;
config.AAD_plot_figure = 1;
config.ISC_select_and_reref = 1;
config.ISC_clean_rawdata = 1;
config.ISC_filt = 1;
config.ISC_forward_model = 1;
config.ISC_same_other = 1;
config.ISC_same_other_chance = 1;
config.ISC_left_right = 1;
config.ISC_plot_figure = 1;
config.ISC_plot_supp_figure = 1;
config.SpE_select_and_reref = 1;
config.SpE_clean_rawdata = 1;
config.SpE_filt = 1;
config.SpE_compute_entropy = 1;
config.SpE_plot_figure = 1;
config.SpE_plot_supp_figure = 1;
config.BTW_plot_AAD_ISC_figure = 1;
config.BTW_plot_AAD_ISC_SpE_supp_figure = 1;

% create data path
DATAPATH = [MAINPATH,'data',filesep];
if ~exist(DATAPATH, 'dir') % create output directory, if necessary
mkdir(DATAPATH)
end


%% Speech Envelope Tracking (AAD)
    %%% select and rereference
        fileout_AAD_select_and_reref = [DATAPATH,'ALLEEG_AAD_select_and_reref'];
        if config.AAD_select_and_reref
            % selects participant specific blocks, rereferences to linked mastoid
            % and removes L04A for symmetry
            [ALLEEG_AAD_select_and_reref] = bjh_01_AAD_01_select_and_reref(MAINPATH,CONFIGPATH);
            save([fileout_AAD_select_and_reref,'.mat'],'ALLEEG_AAD_select_and_reref','-v7.3');
        else
            disp('Loading ALLEEG_AAD_select_and_reref ...');
            load([fileout_AAD_select_and_reref,'.mat'],'ALLEEG_AAD_select_and_reref');
        end

    %%% clean cEEGrid data (clean_rawdata) %%%
        fileout_AAD_clean = [DATAPATH,'ALLEEG_AAD_clean'];
        params.flatline = 'off';
        params.highpass = [0.25, 0.75];
        params.chancrit = 'off';
        params.linenoise = 'off';
        params.burstcrit = 10;
        params.wincrit = 'off';
        params.visualize = 0;
        if config.AAD_clean
            % cleans data using clean_rawdata
            [ALLEEG_AAD_clean] = bjh_01_AAD_02_clean_rawdata(ALLEEG_AAD_select_and_reref,params);
            save([fileout_AAD_clean,'.mat'],'ALLEEG_AAD_clean','-v7.3');
        else
            disp('Loading ALLEEG_AAD_clean ...');
            load([fileout_AAD_clean,'.mat'],'ALLEEG_AAD_clean');
        end

    %%% prepare cEEGrid data for mTRF Toolbox %%%
        fileout_prep_with_asr = [DATAPATH,'ALLEEG_AAD_prep_with_ASR'];
        fileout_prep_without_asr = [DATAPATH,'ALLEEG_AAD_prep_without_ASR'];
        params.lp_filter = 8;
        params.hp_filter = 2;
        params.fs_new = 64;

        if config.AAD_prep 
            % filters, normalizes (only divide by std) and downsamples cEEGrid data
            % epochs data according to blocks [0,600] and concatenates them
                % with ASR
                [ALLEEG_AAD_prep_with_ASR] = bjh_01_AAD_03_filter_normalize_downsample(...
                    ALLEEG_AAD_clean,params);
                save([fileout_prep_with_asr,'.mat'],'ALLEEG_AAD_prep_with_ASR','-v7.3');
                % without ASR
                [ALLEEG_AAD_prep_without_ASR] = bjh_01_AAD_03_filter_normalize_downsample(...
                    ALLEEG_AAD_select_and_reref,params);
                save([fileout_prep_without_asr,'.mat'],'ALLEEG_AAD_prep_without_ASR','-v7.3');
        else
            disp('Loading ALLEEG_AAD_prep_with_ASR ...');
            load([fileout_prep_with_asr,'.mat'],'ALLEEG_AAD_prep_with_ASR');
            disp('Loading ALLEEG_AAD_prep_without_ASR ...');
            load([fileout_prep_without_asr,'.mat'],'ALLEEG_AAD_prep_without_ASR');
        end
    
    
    %%% extract speech envelope %%%
        file_in_audio = dir([MAINPATH,'bjh_cEEGrid_attention',filesep,'stimuli',...
            filesep,'Stream*_fs32k.wav']);
        params.visualize_envelope = 0;

        if config.AAD_envelope 
            % extract speech envelope (absolute hilbert transform, lowpass
            % filter, normalize, downsample)
            bjh_01_AAD_04_extract_envelope_mirkovic2016(file_in_audio,DATAPATH,params);
        end
    
    
    %%% concatenate subject-specifically selected blocks of envelope
        file_in_conc_env = dir([DATAPATH,'Envelope_Mirkovic*.mat']);
        fileout_conc_env_with_ASR = [DATAPATH,'ALLEEG_AAD_with_env_with_ASR'];
        fileout_conc_env_without_ASR = [DATAPATH,'ALLEEG_AAD_with_env_without_ASR'];

        if config.AAD_select_envelope
            % concatenates the subject-specifically selected blocks of the speech
            % envelope and assigns the respective stream as attended or ignored
            % attach envelopes as subfield to ALLEEG structure
                % with ASR
                ALLEEG_AAD_env_with_ASR = bjh_01_AAD_05_select_envelopes(...
                    file_in_conc_env,ALLEEG_AAD_prep_with_ASR);
                save([fileout_conc_env_with_ASR,'.mat'],'ALLEEG_AAD_env_with_ASR','-v7.3');
                % without ASR
                ALLEEG_AAD_env_without_ASR = bjh_01_AAD_05_select_envelopes(...
                    file_in_conc_env,ALLEEG_AAD_prep_without_ASR);
                save([fileout_conc_env_without_ASR,'.mat'],'ALLEEG_AAD_env_without_ASR','-v7.3');
        else
            disp('Loading ALLEEG_AAD_env_with_ASR ...');
            load([fileout_conc_env_with_ASR,'.mat'],'ALLEEG_AAD_env_with_ASR');
            disp('Loading ALLEEG_AAD_env_without_ASR ...');
            load([fileout_conc_env_without_ASR,'.mat'],'ALLEEG_AAD_env_without_ASR');
        end
        

    %%% mTRF Toolbox %%% 
        params.trial_len = 60; % in s
        % percentage of data for independent testing
        params.perc_test = 0.33;
        % number of cross-validation for train/validate and test
        params.nfold_test = 50;
        fileout_val_test_mix_with_ASR = [DATAPATH,...
            'AAD_results_val_test_mix_with_ASR_',num2str(params.trial_len),'s'];
        fileout_val_test_mix_without_ASR = [DATAPATH,...
            'AAD_results_val_test_mix_without_ASR_',num2str(params.trial_len),'s'];
        fileout_val_test_sep_with_ASR = [DATAPATH,...
            'AAD_results_val_test_sep_with_ASR_',num2str(params.trial_len),'s_',...
            num2str(params.perc_test*100),'%_',num2str(params.nfold_test),'_folds'];
        % time lag window is 45 ms
        params.tlag_win = 45;
        % time lag windows are shifted 15 ms (neighbouring time lag windows
        % have 30 ms overlap)
        params.tlag_min = -115:15:575;
        % range from 0 to 300 ms is investigated
        params.tlag_max = params.tlag_min + params.tlag_win;
        % lambda (no regularization at first)
        params.lambdas = 10.^(-5:1:5);
        % indicate whether correlations and accuracies should be visualized
        params.vis_tlag_lambda = 1;

        % validate and test separate
        if config.AAD_val_test_sep
            [AAD_results_val_test_sep_with_ASR] = bjh_01_AAD_06_backward_mTRF_val_test_sep(...
                ALLEEG_AAD_env_with_ASR,params);
            save([fileout_val_test_sep_with_ASR,'.mat'],'AAD_results_val_test_sep_with_ASR');
        else
            load([fileout_val_test_sep_with_ASR,'.mat'],'AAD_results_val_test_sep_with_ASR');
        end

        % validate and test mixed
        if config.AAD_val_test_mix
            % computes correlations between reconstructed and attended as well
            % as reconstructed and ignored speech envelope. Reconstructed
            % envelope is generated using a backward model trained on the
            % attended speech envelope
            % with ASR        
            [AAD_results_val_test_mix_with_ASR] = bjh_01_AAD_06_backward_mTRF_val_test_mix(...
                ALLEEG_AAD_env_with_ASR,params);
            save([fileout_val_test_mix_with_ASR,'.mat'],'AAD_results_val_test_mix_with_ASR','-v7.3');

            % without ASR
            [AAD_results_val_test_mix_without_ASR] = bjh_01_AAD_06_backward_mTRF_val_test_mix(...
                ALLEEG_AAD_env_without_ASR,params);
            save([fileout_val_test_mix_without_ASR,'.mat'],'AAD_results_val_test_mix_without_ASR','-v7.3');
        else
            load([fileout_val_test_mix_with_ASR,'.mat'],'AAD_results_val_test_mix_with_ASR');
            load([fileout_val_test_mix_without_ASR,'.mat'],'AAD_results_val_test_mix_without_ASR');
        end
    
    %%% plot speech envelope figure %%%
        fileout_plot_AAD_figure = [DATAPATH,'Figure_2'];

        if config.AAD_plot_figure
            bjh_01_AAD_07_plot_AAD_figure(fileout_plot_AAD_figure,AAD_results_val_test_mix_with_ASR,...
                AAD_results_val_test_mix_without_ASR,AAD_results_val_test_sep_with_ASR,params);
        end

        
%% Intersubject Correlation (ISC)

    %%% import cEEGrid data %%% 
        fileout_ISC_select_and_reref = [DATAPATH,'ALLEEG_ISC_select_and_reref'];

        if config.ISC_select_and_reref
            % epochs according to 10-minute blocks [-5, 605], only keeps first block, 
            % rereferences to linked mastoid, and removes L04A for symmetry
            [ALLEEG_ISC_select_and_reref] = bjh_02_ISC_01_select_and_reref(MAINPATH,CONFIGPATH);
            save([fileout_ISC_select_and_reref,'.mat'],'ALLEEG_ISC_select_and_reref','-v7.3');
        else
            disp('Loading ALLEEG_ISC_select_and_reref ...');
            load([fileout_ISC_select_and_reref,'.mat'],'ALLEEG_ISC_select_and_reref');
        end

    %%% clean cEEGrid data (clean_rawdata) %%%
        fileout_ISC_clean = [DATAPATH,'ALLEEG_ISC_clean'];
        params.flatline = 'off';
        params.highpass = [0.25, 0.75];
        params.chancrit = 'off';
        params.linenoise = 'off';
        params.burstcrit = 10;
        params.wincrit = 'off';
        params.visualize = 0;

        if config.ISC_clean_rawdata
            % cleans data using clean_rawdata
            [ALLEEG_ISC_clean] = bjh_02_ISC_02_clean_rawdata(ALLEEG_ISC_select_and_reref,params);
            save([fileout_ISC_clean,'.mat'],'ALLEEG_ISC_clean','-v7.3');
        else
            disp('Loading ALLEEG_ISC_clean ...');
            load([fileout_ISC_clean,'.mat'],'ALLEEG_ISC_clean');
        end

    %%% filter and downsample data %%%
        file_out_filt_with_asr = [DATAPATH,'ALLEEG_ISC_filt'];
        params.lp_filter = 40;
        params.hp_filter = 1;
        params.fs_new = 250;

        if config.ISC_filt
            % filters and downsamples cEEGrid data
            [ALLEEG_ISC_filt] = bjh_02_ISC_03_filt_downsample(ALLEEG_ISC_clean,params);
            save([file_out_filt_with_asr,'.mat'],'ALLEEG_ISC_filt','-v7.3');
        else
            disp('Loading ALLEEG_filt_with_ASR ...');
            load([file_out_filt_with_asr,'.mat'],'ALLEEG_ISC_filt');
        end
    
    
    %%% ISC Forward Model %%%
        fileout_ISC_fmodel = [DATAPATH,'ISC_fmodel'];
        params.gamma = 0.4;
        params.Ncomp = 3;

        if config.ISC_forward_model
            % computes forward model (condition-independent, left, and right)
            [A_cond_indep,A_left,A_right] = bjh_02_ISC_04_fmodel(ALLEEG_ISC_filt,params);
            save([fileout_ISC_fmodel,'.mat'],'A_cond_indep','A_left','A_right','-v7.3');
        else
            disp('Loading A_cond_indep, A_left, A_right ...');
            load([fileout_ISC_fmodel,'.mat'],'A_cond_indep','A_left','A_right');
        end

    %%% Compute ISC_same / ISC_other %%%
        fileout_ISC_same_other = [DATAPATH,'ISC_results_same_other'];

        if config.ISC_same_other
            % compute condition dependent ISC analysis (ISC_same, ISC_other)
            [ISC_same_ind,ISC_other_ind] = bjh_02_ISC_04_same_other(ALLEEG_ISC_filt,params);
            fprintf('\n');
            save([fileout_ISC_same_other,'.mat'],'ISC_same_ind','ISC_other_ind','-v7.3');
        else
            disp('Loading ISC_same_ind, ISC_other_ind ...');
            load([fileout_ISC_same_other,'.mat'],'ISC_same_ind','ISC_other_ind');  
        end
    
        %%% Compute Chance Level - ISC_same / ISC_other
            fileout_ISC_same_other_chance = [DATAPATH,'ISC_results_same_other_chance'];
            imax = 100;
            ISC_same_ind_chance = zeros(size(ALLEEG_ISC_filt(1).data,1),...
                size(ALLEEG_ISC_filt,2),imax);
            ISC_other_ind_chance = ISC_same_ind_chance;

            if config.ISC_same_other_chance
                for i = 1:imax
                    fprintf('Iteration %d: ',i);
                    % compute condition dependent ISC analysis (ISC_same, ISC_other)
                    [ISC_same_ind_chance(:,:,i),ISC_other_ind_chance(:,:,i)] = ...
                        bjh_02_ISC_04_same_other(ALLEEG_ISC_filt,params,1);
                    fprintf('\n');
                end
                save([fileout_ISC_same_other_chance,'.mat'],...
                    'ISC_same_ind_chance','ISC_other_ind_chance','-v7.3');
            else
                load([fileout_ISC_same_other_chance,'.mat'],...
                    'ISC_same_ind_chance','ISC_other_ind_chance');
                disp('Loading ISC_same_ind_chance, ISC_other_ind_chance ...');
            end
    
    %%% Compute ISC_left and ISC_right %%%
        fileout_ISC_left_right = [DATAPATH,'ISC_results_left_right'];

        if config.ISC_left_right
            % compute ISC left and ISC right
            [ISC_left_ind,ISC_right_ind,AUC_left,AUC_right,AUC_95_prctile] = ...
                bjh_02_ISC_04_left_right(ALLEEG_ISC_filt,params);
            save([fileout_ISC_left_right,'.mat'],...
                'ISC_left_ind','ISC_right_ind','AUC_left',...
                'AUC_right','AUC_95_prctile','-v7.3');
        else
            disp('Loading ISC_left_ind, ISC_right_ind, AUC_left, AUC_right, AUC_95_prctile ...');
            load([fileout_ISC_left_right,'.mat'],...
                'ISC_left_ind','ISC_right_ind','AUC_left',...
                'AUC_right','AUC_95_prctile');
        end
        
    %%% Plot ISC Figure %%%
        fileout_plot_ISC_figure = [DATAPATH,'Figure_3'];

        if config.ISC_plot_figure
            bjh_02_ISC_05_plot_ISC_figure(CONFIGPATH,fileout_plot_ISC_figure,ISC_same_ind,ISC_other_ind,ALLEEG_ISC_filt,...
                ISC_same_ind_chance,ISC_other_ind_chance,ISC_left_ind,ISC_right_ind,A_cond_indep,params);
        end

    %%% Plot ISC Supplementary Figure %%%
        fileout_plot_ISC_supp_figure = [DATAPATH,'Supp_Figure_1'];

        if config.ISC_plot_supp_figure
            bjh_02_ISC_05_plot_ISC_supp_fig(fileout_plot_ISC_supp_figure,ALLEEG_ISC_filt,...
                A_left,A_right,params)
        end
        
%% Spectral Entropy (SpE)

    %%% select blocks and rereference cEEGrid data %%%
        fileout_SpE_select_and_reref = [DATAPATH,'ALLEEG_SpE_select_and_reref'];

        if config.SpE_select_and_reref
            % epochs according to 10-minute blocks [-5, 605], only keeps pre-
            % selected blocks and concatenates them to get continuous data,
            % rereferences to linked mastoid, and removes L04A for symmetry
            [ALLEEG_SpE_select_and_reref] = bjh_03_SpE_01_select_and_reref(MAINPATH,CONFIGPATH);
            save([fileout_SpE_select_and_reref,'.mat'],'ALLEEG_SpE_select_and_reref','-v7.3');
        else
            disp('Loading ALLEEG_SpE_select_and_reref ...');
            load([fileout_SpE_select_and_reref,'.mat'],'ALLEEG_SpE_select_and_reref');
        end


    %%% clean cEEGrid data (clean_rawdata)
        fileout_SpE_clean = [DATAPATH,'ALLEEG_SpE_clean'];
        params.flatline = 'off';
        params.highpass = [0.25, 0.75];
        params.chancrit = 0.6;
        params.linenoise = 'off';
        params.burstcrit = 10;
        params.wincrit = 'off';
        params.visualize = 0;

        if config.SpE_clean_rawdata
            % cleans data using clean_rawdata
            [ALLEEG_SpE_clean] = bjh_03_SpE_02_clean_rawdata(ALLEEG_SpE_select_and_reref,params);
            save([fileout_SpE_clean,'.mat'],'ALLEEG_SpE_clean','-v7.3');
        else
            disp('Loading ALLEEG_SpE_clean ...');
            load([fileout_SpE_clean,'.mat'],'ALLEEG_SpE_clean');
        end


    %%% filter cEEGrid data %%%
        fileout_SpE_filt = [DATAPATH,'ALLEEG_SpE_filt'];
        params.lp_filter = 40;
        params.hp_filter = 1;

        if config.SpE_filt
            % filters and epochs data according to blocks [0,600] and concatenates them
            [ALLEEG_SpE_filt] = bjh_03_SpE_03_filt(ALLEEG_SpE_clean,params);
            save([fileout_SpE_filt,'.mat'],'ALLEEG_SpE_filt','-v7.3');
        else
            disp('Loading ALLEEG_SpE_filt ...');
            load([fileout_SpE_filt,'.mat'],'ALLEEG_SpE_filt');
        end
    
    %%% Compute spectral entropy %%%
        fileout_SpE_entropy = [DATAPATH,'SpE_results'];
        params.seg_length = 60; % in s
        params.entr_freq_min = 8; % in Hz
        params.entr_freq_max = 32; % in Hz

        if config.SpE_compute_entropy
            [SpE_results] = bjh_03_SpE_04_compute_entropy(ALLEEG_SpE_filt,params);
            save([fileout_SpE_entropy,'.mat'],'SpE_results','-v7.3');
        else
            disp('Loading SpE_results ...');
            load([fileout_SpE_entropy,'.mat'],'SpE_results')
        end
        
    %%% Plot Spectral Entropy Figure %%%
        fileout_SpE_plot_figure = [DATAPATH,'Figure_4'];

        if config.SpE_plot_figure
            bjh_03_SpE_05_plot_entropy_figure(fileout_SpE_plot_figure,SpE_results,params)
        end

    %%% Plot Spectral Entropy Supplementary Figure %%%
        fileout_SpE_plot_supp_figure = [DATAPATH,'Supp_Figure_2'];

        if config.SpE_plot_supp_figure
            bjh_03_SpE_05_plot_entropy_supp_figure(fileout_SpE_plot_supp_figure,SpE_results)
        end
        
        
%% Between Measures
    
    %%% plot AAD ~ ISC figure %%%
        fileout_BTW_plot_AAD_ISC_figure = [DATAPATH,'Figure_5'];

        if config.BTW_plot_AAD_ISC_figure == 1
            bjh_04_BTW_01_plot_AAD_ISC_figure(fileout_BTW_plot_AAD_ISC_figure,...
                ISC_same_ind,ISC_other_ind,AAD_results_val_test_mix_with_ASR,params);
        end
    
    % plot AAD ~ SpE and ISC ~ SpE supplementary figure
        fileout_plot_AAD_ISC_SpE_supp_figure = [DATAPATH,'Supp_Figure_3'];

        if config.BTW_plot_AAD_ISC_SpE_supp_figure == 1
            bjh_04_BTW_01_plot_AAD_ISC_SpE_supp_figure(fileout_plot_AAD_ISC_SpE_supp_figure,...
                AAD_results_val_test_mix_with_ASR,SpE_results,ISC_same_ind,ISC_other_ind,params)  
        end
    
    
    

