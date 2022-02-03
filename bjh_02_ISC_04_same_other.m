function [ISC_same_ind,ISC_other_ind] = bjh_02_ISC_04_same_other(ALLEEG_ISC_filt,params,chance)
%% computes condition dependent ISC (ISC_same, ISC_other)
% input:    ALLEEG_ISC_filt:  structure array containing EEG structures from each
%                       partcipant (similar to ALLEEG in eeglab)
%           params:     structure reguöarzation parameter gamma
%           chance:     boolean whether chance level should be computed
% 
% output:   ISC_same_ind:    [double] array (16 ISC x 36 subjects) ISC components
%                               computed with subjects attending to the same story 
%           ISC_other_ind:   [double] array (16 ISC x 36 subjects) ISC components
%                               computed with subjects attending to the other story 
%
% author: 
%   Original code: (c) Lucas C Parra, parra@ccny.cuny.edu (isceeg.m
%                   (https://www.parralab.org/isc/))
%   Adaptation: Björn Holtze
% date: 14.01.2022

if nargin == 2
    chance = 0;
end

ISC_same_ind = zeros(size(ALLEEG_ISC_filt(1).data,1),size(ALLEEG_ISC_filt,2));
ISC_other_ind = zeros(size(ALLEEG_ISC_filt(1).data,1),size(ALLEEG_ISC_filt,2));

% preallocate space to data array X (samples x channels x subjects)
X = zeros(size(ALLEEG_ISC_filt(1).data,2),size(ALLEEG_ISC_filt(1).data,1),...
    size(ALLEEG_ISC_filt,2));

% create matrix X (samples x channels x subjects)
for i = 1:size(ALLEEG_ISC_filt,2)
    if chance == 0
        X(:,:,i) = permute(ALLEEG_ISC_filt(i).data,[2,1]);
    elseif chance == 1
        X(:,:,i) = circshift(permute(ALLEEG_ISC_filt(i).data,[2,1]),...
            randi(size(permute(ALLEEG_ISC_filt(i).data,[2,1]),1)));
    end
end

%% Compute cross-covariance for all participants
    % T samples, D channels, N subjects
    [~,D,N] = size(X);

    % compute cross-covariance between all subjects i and j
    Rij = permute(reshape(cov(X(:,:)),[D,N,D,N]),[1,3,2,4]);
    
for subj = 1:size(ALLEEG_ISC_filt,2)
    fprintf('%d ', subj);
    
    %% compute projection vector W (without to-be-correlated participant)
        % exclude to-be-correlated subject from projection vector calculation
        X_wo_subj = X(:,:,1:size(ALLEEG_ISC_filt,2)~= subj);

        % T samples, D channels, N subjects
        [~,D_wo_subj,N_wo_subj] = size(X_wo_subj);

        % compute cross-covariance between all subjects i and j
        Rij_wo_subj = permute(reshape(cov(X_wo_subj(:,:)),[D_wo_subj,N_wo_subj,...
            D_wo_subj,N_wo_subj]),[1,3,2,4]);

        % compute within- and between-subject covariances
        Rw_wo_subj =               1/N_wo_subj * ...
            sum(Rij_wo_subj(:,:,1:N_wo_subj+1:N_wo_subj*N_wo_subj),3); % pooled over all subjects
        Rb_wo_subj = 1/(N_wo_subj-1)/N_wo_subj  * ...
            (sum(Rij_wo_subj(:,:,:),3) - N_wo_subj*Rw_wo_subj); % pooled over all pairs of subjects

        % shrinkage regularization of Rw
        Rw_reg_wo_subj = (1-params.gamma)*Rw_wo_subj + ...
            params.gamma*mean(eig(Rw_wo_subj))*eye(size(Rw_wo_subj));

        % compute correlated components W using regularized Rw, sort components by ISC
        [W_wo_subj,ISC_wo_subj] = eig(Rb_wo_subj,Rw_reg_wo_subj);
        [~,idx] = sort(diag(ISC_wo_subj),'descend');
        W_wo_subj = W_wo_subj(:,idx);
    
        
    %% ISC same
        % Compute ISC resolved by subject, see Cohen et al.
        idx_same = [ALLEEG_ISC_filt.attended_ch] == ALLEEG_ISC_filt(subj).attended_ch;
        N_same = sum(idx_same);
        
        % compute within-subject covariance (same attended side)
        Rw_same = 0;
        for j=1:N
            % if subject j attended to the same story as subj
            if idx_same(j) == 1
                if subj~=j
                    Rw_same = Rw_same+1/(N_same-1)*(Rij(:,:,subj,subj)+Rij(:,:,j,j));
                end
            end
        end
        % compute within-subject covariance (same attended side)
        Rb_same = 0;
        for j=1:N
            % if subject j attended to the same story as subj
            if idx_same(j) == 1
                if subj~=j
                    Rb_same = Rb_same+1/(N_same-1)*(Rij(:,:,subj,j)+Rij(:,:,j,subj));
                end
            end
        end
        % compute ISC same
        ISC_same_ind(:,subj) = diag(W_wo_subj'*Rb_same*W_wo_subj)./diag(W_wo_subj'*Rw_same*W_wo_subj);
    
     %% ISC other
        % Compute ISC resolved by subject, see Cohen et al.
        idx_other = [ALLEEG_ISC_filt.attended_ch] ~= ALLEEG_ISC_filt(subj).attended_ch;
        N_other = sum(idx_other);
        
        % compute within-subject covariance (same attended side)
        Rw_other = 0;
        for j=1:N
            % if subject j attended to the other story than subj
            if idx_other(j) == 1
                if subj~=j
                    Rw_other = Rw_other+1/(N_other-1)*(Rij(:,:,subj,subj)+Rij(:,:,j,j));
                end
            end
        end
        % compute within-subject covariance (same attended side)
        Rb_other = 0;
        for j=1:N
            % if subject j attended to the same story as subj
            if idx_other(j) == 1
                if subj~=j
                    Rb_other = Rb_other+1/(N_other-1)*(Rij(:,:,subj,j)+Rij(:,:,j,subj));
                end
            end
        end
        % compute ISC same
        ISC_other_ind(:,subj) = diag(W_wo_subj'*Rb_other*W_wo_subj)./diag(W_wo_subj'*Rw_other*W_wo_subj);

end
end

