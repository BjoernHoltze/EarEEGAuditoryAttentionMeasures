function [ISC_left_ind,ISC_right_ind,AUC_left,AUC_right,AUC_95_prctile] = ...
    bjh_02_ISC_04_left_right(ALLEEG_ISC_filt,params)
%% computes ISC left and ISC right
% input:    ALLEEG_filt_with_ASR:  structure array containing EEG structures from each
%                       partcipant (similar to ALLEEG in eeglab)
%           params:     structure reguöarzation parameter gamma
% 
% output:   ISC_left_ind:   [double] array (16 ISC x 36 subjects) ISC components
%                           computed with subjects attending to the left story 
%           ISC_right_ind:  [double] array (16 ISC x 36 subjects) ISC components
%                           computed with subjects attending to the right story 
%           AUC_left:       Area under the receiver operator curve for ISC_left
%           AUC_right:      Area under the receiver operator curve for ISC_right
%           AUC_95_prctile: Area under the receiver operator curve for
%                           randomly shuffed labels (chance level)
%
% author: 
%   Original code: (c) Lucas C Parra, parra@ccny.cuny.edu (isceeg.m
%                   (https://www.parralab.org/isc/))
%   Adaptation: Björn Holtze
% date: 14.01.2022

ISC_left_ind = zeros(size(ALLEEG_ISC_filt(1).data,1),size(ALLEEG_ISC_filt,2));
ISC_right_ind = zeros(size(ALLEEG_ISC_filt(1).data,1),size(ALLEEG_ISC_filt,2));
attended_ch = [ALLEEG_ISC_filt.attended_ch];

% preallocate space to data array X (samples x channels x subjects)
X = zeros(size(ALLEEG_ISC_filt(1).data,2),size(ALLEEG_ISC_filt(1).data,1),...
    size(ALLEEG_ISC_filt,2));

% create matrix X (samples x channels x subjects)
for i = 1:size(ALLEEG_ISC_filt,2)
    X(:,:,i) = permute(ALLEEG_ISC_filt(i).data,[2,1]);
end

% Compute cross-covariance (D channels, N subjects)
[~,D,N] = size(X);

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D,N,D,N]),[1,3,2,4]);

for subj = 1:size(ALLEEG_ISC_filt,2)
    fprintf('%d ', subj);
    
    %% remove to be correlated participant from X
    X_wo_subj = X(:,:,1:size(ALLEEG_ISC_filt,2)~= subj);
        
    %% ISC left 
    
        %% Projection Vector (only considering those attending to the left story)
            % extract data only from those attending to the left
            % (to-be-correlated participants is excluded already)
            X_left = X_wo_subj(:,:,attended_ch(1:size(ALLEEG_ISC_filt,2) ~= subj) == 1);
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
        
        %% compute individual ISC left value of the to-be-correlated participant
        % Compute ISC resolved by subject, see Cohen et al.
        idx_left = attended_ch == 1;
        
        % compute within-subject covariance (attended left)
        Rw_left_ind = 0;
        for j=1:N
            % if subject j attended to the left story
            if idx_left(j) == 1
                if subj~=j
                    Rw_left_ind = Rw_left_ind+1/(N_left-1)*(Rij(:,:,subj,subj)+Rij(:,:,j,j));
                end
            end
        end
        % compute within-subject covariance (attended left)
        Rb_left_ind = 0;
        for j=1:N
            % if subject j attended to the left story
            if idx_left(j) == 1
                if subj~=j
                    Rb_left_ind = Rb_left_ind+1/(N_left-1)*(Rij(:,:,subj,j)+Rij(:,:,j,subj));
                end
            end
        end
        % compute ISC left ind
        ISC_left_ind(:,subj) = diag(W_left'*Rb_left_ind*W_left)./diag(W_left'*Rw_left_ind*W_left);
    
    %% ISC right
    
        %% Projection Vector (only considering those attending to the right story)
            % extract data only from those attending to the right
            % (to-be-correlated participants is excluded already)
            X_right = X_wo_subj(:,:,attended_ch(1:size(ALLEEG_ISC_filt,2) ~= subj) == 2);
            % T samples, D channels, N subjects
            [~,D_right,N_right] = size(X_right);

            % compute cross-covariance between all subjects attending to the right story
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
        
        %% compute individual ISC left value of the to-be-correlated participant
        % Compute ISC resolved by subject, see Cohen et al.
        idx_right = attended_ch == 2;
        
        % compute within-subject covariance (attended left)
        Rw_right_ind = 0;
        for j=1:N
            % if subject j attended to the left story
            if idx_right(j) == 1
                if subj~=j
                    Rw_right_ind = Rw_right_ind+1/(N_right-1)*(Rij(:,:,subj,subj)+Rij(:,:,j,j));
                end
            end
        end
        % compute within-subject covariance (attended left)
        Rb_right_ind = 0;
        for j=1:N
            % if subject j attended to the left story
            if idx_right(j) == 1
                if subj~=j
                    Rb_right_ind = Rb_right_ind+1/(N_right-1)*(Rij(:,:,subj,j)+Rij(:,:,j,subj));
                end
            end
        end
        % compute ISC left ind
        ISC_right_ind(:,subj) = diag(W_right'*Rb_right_ind*W_right)./diag(W_right'*Rw_right_ind*W_right);
    

end
fprintf('\n ');

%% Classification (AUC of ROC)

% calculate ISC_left_score and ISC_right_score
ISC_left_ind_score = sum(ISC_left_ind(1:params.Ncomp,:),1);
ISC_right_ind_score = sum(ISC_right_ind(1:params.Ncomp,:),1);

[~,~,~,AUC_left] = perfcurve(attended_ch,ISC_left_ind_score,1);
[~,~,~,AUC_right] = perfcurve(attended_ch,ISC_right_ind_score,2);

for i = 1:1000
    attended_ch_chance = attended_ch(randperm(size(attended_ch,2)));
    [~,~,~,AUC_left_chance(i)] = perfcurve(attended_ch_chance,ISC_left_ind_score,1);
    [~,~,~,AUC_right_chance(i)] = perfcurve(attended_ch_chance,ISC_right_ind_score,2);
end

    AUC_95_prctile = prctile([AUC_left_chance,AUC_right_chance],95);

end

