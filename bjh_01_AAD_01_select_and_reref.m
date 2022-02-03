function [ALLEEG_select_and_reref] = bjh_01_AAD_01_select_and_reref(MAINPATH,CONFIGPATH)
%% selects relevant blocks, and references to linked mastoids
% epochs according to 10 minute blocks ([-10, 610]), removes DC
% offset (first 500 ms of each epoch), keeps only the previously selected 
% 10 minute blocks (3 blocks per participants) and converts epochs data
% back to continuous data. Lastly, the data is rereferenced to linked 
% mastoid and channel L04A is removed to have same number of channels
% on the left and right side.
% 
% input:    MAINPATH:   [string] path used to load .set files
%           CONFIGPATH: [string] path where subj_info.mat is stored
% 
% output:   ALLEEG_imp: structure array containing EEG structures as 
%                       subfields
%
% author: Björn Holtze
% date: 14.01.22

% load subject information (id, attended channel, selected blocks)
load([CONFIGPATH,'subj_info.mat'],'subj_info');

for s = 1:size(subj_info,1)
    disp(['Processing subject ',num2str(s),' ...']);
    
    % load dataset
    if length(num2str(s)) == 1
        subj = ['00',num2str(s)];
    elseif length(num2str(s)) == 2
        subj = ['0',num2str(s)];
    end
    EEG = pop_loadset('filename',['sub-',subj,'_task-AttendedSpeakerParadigmcEEGridAttention_eeg.set'],...
        'filepath',[MAINPATH,'bjh_cEEGrid_attention',filesep,'sub-',subj,filesep,'eeg',filesep]);   
    
    % epoch to extract 10 minute blocks
    EEG = pop_epoch(EEG,{'StartTrigger'},[-5  605],'epochinfo','yes');
    EEG.epochs_before_selection = EEG.trials;
    
    % remove DC offset from  each epoch (to make highpass filter faster) 
    EEG = pop_rmbase(EEG,[-5000, -4500]);
    
    % select relevant epochs (different for each participant)
    EEG = pop_selectevent(EEG,'type',{'StartTrigger'},'epoch',subj_info(s).selected_bl,...
        'deleteevents','off','deleteepochs','on','invertepochs','off');
    
    % concatenate epoched to continuous data
    EEG = eeg_epoch2continuous(EEG);
    
    % rereference to linked mastoid
    refs = find(strcmpi({EEG.chanlocs.labels},'L04B'));
    EEG.data(refs,:) = 0.5*EEG.data(refs,:);
    EEG = pop_reref(EEG,refs);
    EEG.ref = 'linked mastoid';

    % remove L04B for symmetry
    EEG = pop_select(EEG,'nochannel',{'L04A'});
    
    % add subject subj_info
    EEG.setname = subj_info(s).subj_id;
    EEG.attended_ch = subj_info(s).attended_ch;
    EEG.selected_bl = subj_info(s).selected_bl;
    
    % store EEG in ALLEEG(s)
    ALLEEG_select_and_reref(s) = EEG;
    clear EEG;
   
end
end

