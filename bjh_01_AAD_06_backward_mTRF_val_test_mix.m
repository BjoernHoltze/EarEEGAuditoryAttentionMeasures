function [results] = bjh_01_AAD_06_backward_mTRF_val_test_mix(ALLEEG_env,params)
%% computes the correlation as a function of time lag windows and lambdas
% In particular the correlation between the attended and the reconstructed
% speech envelope as well as the correlation between the ignored and the
% reconstructed speech envelope is computed for different time lag windows.
% The reconstructed envelope is generated using a backward model trained on
% the attended speech envelope in a cross-validation procedure (mTRFattncrossval.m)
% For that, the continuous data of a participant is first partitioned into trials 
% using the mTRFpartition.m function. This procedure is done for each participant 
% separately. 
% 
% input:    ALLEEG_env:     structure containing preprocessed EEG data and
%                           speech envelope as subfiled (similar to ALLEEG
%                           in EEGLAB)
%           file_out_time_lag: filename including full path of mat files to
%                           be stored (without .mat extension)
%           params:         structure containing the following relevant
%                           subfield: 
%                           - trial_len: length of a trial in s
%                           - tlag_min: vector containing lower boundaries 
%                           of time lag window
%                           - tlag_max: vector containing upper boundaries 
%                           of time lag window
%                           - lambdas: scalar specifying the lambda
%                           used for regularization
% 
% output:   acc:            array containing classification accuracies
%                           (subject x number of time lag windows x lambdas)
%           corr:           structure containing subfields "att" and "ign"
%                           each containing matrix of correlations with 
%                           respective envelope
%                           (subject x number of time lag windows x lambdas x trials)     
%           err:            structure containing subfields "att" and "ign"
%                           each containing matrix with prediction error of 
%                           respective envelope
%                           (subject x number of time lag windows x lambdas x trials)
% 
% author: Björn Holtze 
% date: 22.04.21

    nb_trials = size(ALLEEG_env(1).data,2)/(ALLEEG_env(1).srate*params.trial_len);
    results = struct;
    results.acc.cv = zeros(size(ALLEEG_env,2),length(params.tlag_min),length(params.lambdas));
    results.acc.ind = zeros(size(ALLEEG_env,2),1);
    results.acc.ga = zeros(size(ALLEEG_env,2),1);
    results.corr.att = zeros(size(ALLEEG_env,2),length(params.tlag_min),length(params.lambdas),nb_trials);
    results.corr.ign = zeros(size(ALLEEG_env,2),length(params.tlag_min),length(params.lambdas),nb_trials);
    results.err.att = zeros(size(ALLEEG_env,2),length(params.tlag_min),length(params.lambdas),nb_trials);
    results.err.ign = zeros(size(ALLEEG_env,2),length(params.tlag_min),length(params.lambdas),nb_trials);

for s = 1:size(ALLEEG_env,2)
    disp(['Processing subject ',num2str(s),' ...']);
    
    % partition cEEGrid data and envelope
        % attended 
        [stim_att_train{s,1},resp_train{s,1}] = mTRFpartition(ALLEEG_env(s).env_att',...
            double(ALLEEG_env(s).data)',nb_trials);
        % ignored
        [stim_ign_train{s,1},~] = mTRFpartition(ALLEEG_env(s).env_ign',...
            double(ALLEEG_env(s).data)',nb_trials);
    
    % perform cross-validation approach for different time lag windows
    for t = 1:size(params.tlag_min,2)
        % attended stream
        
        [~,cv_att,cv_ign] = mTRFattncrossval(stim_att_train{s,1},...
            stim_ign_train{s,1},resp_train{s,1},ALLEEG_env(s).srate,...
            -1,params.tlag_min(t),params.tlag_max(t),params.lambdas,...
            'method','Tikhonov','type','multi','corr','Pearson',...
            'error','mse','zeropad',1,'fast',1,'verbose',0);
        
        % calculate classification accuracy as the number trials where 
        % correlation of attended is larger then correlation of ignored
        % divided by the number of total trials
        results.acc.cv(s,t,:) = (sum(cv_att.r > cv_ign.r)/nb_trials)*100;
        results.corr.att(s,t,:,:) = cv_att.r';
        results.corr.ign(s,t,:,:) = cv_ign.r';
        results.err.att(s,t,:,:) = cv_att.err';
        results.err.ign(s,t,:,:) = cv_ign.err';
    end
    
end

%% Identify Optimal Parameter
    % Individually chosen
    for s = 1:size(ALLEEG_env,2)
        acc_cv_ind = squeeze(results.acc.cv(s,:,:));
        % determine maximum decoding accuracy
        max_acc_val = max(acc_cv_ind(:));
        err_ind = squeeze(mean(results.err.att(s,:,:,:),4));
        % put error to 1 wherever decoding accuracy is not at maximum
        err_ind(acc_cv_ind < max_acc_val) = 1;
        % get index where error is minimal only considering the values
        % where decoding accuracy is maximal
        [min_err_row,min_err_row_idx] = min(err_ind);
        [~,min_err_lambda_idx] = min(min_err_row);
        min_err_lag_idx = min_err_row_idx(min_err_lambda_idx);
        % store optimal time lag window, lambda and the corresponding accuracy
        results.optimal_params.ind_tlag(s,:) = [params.tlag_min(min_err_lag_idx),...
            params.tlag_max(min_err_lag_idx)];
        results.optimal_params.ind_lambda(s,:) = params.lambdas(min_err_lambda_idx);
        results.acc.ind(s,:) = max_acc_val;
    end

	% mutually chosen
    [max_row_ga, idx_row_ga] = max(squeeze(mean(results.acc.cv,1)));
    [~, max_acc_idx_lambda_ga] = max(max_row_ga);
    max_acc_idx_tlag_ga = idx_row_ga(max_acc_idx_lambda_ga);
    results.optimal_params.ga_tlag = [params.tlag_min(max_acc_idx_tlag_ga),...
        params.tlag_max(max_acc_idx_tlag_ga)];
    results.optimal_params.ga_lambda = params.lambdas(max_acc_idx_lambda_ga);
    
    for s = 1:size(ALLEEG_env,2)
        results.acc.ga(s,:) = results.acc.cv(s,max_acc_idx_tlag_ga,max_acc_idx_lambda_ga);
    end
end


