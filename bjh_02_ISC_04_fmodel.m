function [A_cond_indep,A_left,A_right] = bjh_02_ISC_04_fmodel(ALLEEG_ISC_filt,params)
%% computes condition independent ISC (to obtain forward model)
% input:    ALLEEG_filt:structure array containing EEG structures from each
%                       partcipant (similar to ALLEEG in eeglab)
%           params:     regularzation parameter gamma and number
%                       of components to be summed for ISC score
% 
% output:   A_cond_indep:   Condition independent ISC forward model
%           A_left:         ISC forward model (left)
%           A_right:        ISC forward model (right)
%
% author: 
%   Original code: (c) Lucas C Parra, parra@ccny.cuny.edu (isceeg.m
%                   (https://www.parralab.org/isc/))
% Adaptation: Björn Holtze
% date: 14.01.2022

%% ISC Forward Model (over all particpants)
    X = zeros(size(ALLEEG_ISC_filt(1).data,2),size(ALLEEG_ISC_filt(1).data,1),...
        size(ALLEEG_ISC_filt,2));
    for s = 1:size(ALLEEG_ISC_filt,2)
        X(:,:,s) = permute(ALLEEG_ISC_filt(s).data,[2,1]);
    end

    % T samples, D channels, N subjects
    [~,D,N] = size(X);  

    % compute cross-covariance between all subjects i and j
    Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

    % compute within- and between-subject covariances
    Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
    Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

    % shrinkage regularization of Rw
    Rw_reg = (1-params.gamma)*Rw + params.gamma*mean(eig(Rw))*eye(size(Rw));

    % compute correlated components W using regularized Rw, sort components by ISC
    [W_cond_indep,ISC_cond_indep] = eig(Rb,Rw_reg); 
    [~,indx] = sort(diag(ISC_cond_indep),'descend');
    W_cond_indep = W_cond_indep(:,indx);

    % compute forward model ("scalp projections") A
    A_cond_indep = Rw*W_cond_indep/(W_cond_indep'*Rw*W_cond_indep);

    
%% ISC Forward Model (ISC_left)
    attended_ch = [ALLEEG_ISC_filt.attended_ch];
    X_left = X(:,:,attended_ch(1:size(ALLEEG_ISC_filt,2)) == 1);
    % T samples, D channels, N subjects
    [~,D_left,N_left] = size(X_left);

    % compute cross-covariance between all subjects i and j
    Rij_left = permute(reshape(cov(X_left(:,:)),[D_left,N_left,...
        D_left,N_left]),[1,3,2,4]);

    % compute within- and between-subject covariances
    Rw_left =               1/N_left * ...
        sum(Rij_left(:,:,1:N_left+1:N_left*N_left),3); % pooled over all subjects
    Rb_left = 1/(N_left-1)/N_left  * ...
        (sum(Rij_left(:,:,:),3) - N_left*Rw_left); % pooled over all pairs of subjects

    % shrinkage regularization of Rw
    Rw_reg_left = (1-params.gamma)*Rw_left + ...
        params.gamma*mean(eig(Rw_left))*eye(size(Rw_left));

    % compute correlated components W using regularized Rw
    [W_left,ISC_left] = eig(Rb_left,Rw_reg_left);
    [~,left_idx] = sort(diag(ISC_left),'descend');
    W_left = W_left(:,left_idx);
    
    % compute forward model ("scalp projections") A
    A_left = Rw_left*W_left/(W_left'*Rw_left*W_left);
    
    
%% ISC Forward Model (ISC_right)
    attended_ch = [ALLEEG_ISC_filt.attended_ch];
    X_right = X(:,:,attended_ch(1:size(ALLEEG_ISC_filt,2)) == 2);
    % T samples, D channels, N subjects
    [~,D_right,N_right] = size(X_right);

    % compute cross-covariance between all subjects i and j
    Rij_right = permute(reshape(cov(X_right(:,:)),[D_right,N_right,...
        D_right,N_right]),[1,3,2,4]);

    % compute within- and between-subject covariances
    Rw_right =               1/N_right * ...
        sum(Rij_right(:,:,1:N_right+1:N_right*N_right),3); % pooled over all subjects
    Rb_right = 1/(N_right-1)/N_right  * ...
        (sum(Rij_right(:,:,:),3) - N_right*Rw_right); % pooled over all pairs of subjects

    % shrinkage regularization of Rw
    Rw_reg_right = (1-params.gamma)*Rw_right + ...
        params.gamma*mean(eig(Rw_right))*eye(size(Rw_right));

    % compute correlated components W using regularized Rw
    [W_right,ISC_right] = eig(Rb_right,Rw_reg_right);
    [~,right_idx] = sort(diag(ISC_right),'descend');
    W_right = W_right(:,right_idx);
    
    % compute forward model ("scalp projections") A
    A_right = Rw_right*W_right/(W_right'*Rw_right*W_right);


end

