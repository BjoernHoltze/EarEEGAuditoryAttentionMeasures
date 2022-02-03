function bjh_04_BTW_01_plot_AAD_ISC_SpE_supp_figure(fileout_plot_AAD_ISC_SpE_supp_figure,...
            AAD_results_val_test_mix_with_ASR,SpE_results,ISC_same_ind,ISC_other_ind,params)   
%% Computation

    %% corr_att - corr_ign
    corr_att = zeros(36,30);
    corr_ign = zeros(36,30);
    for s = 1:size(AAD_results_val_test_mix_with_ASR.corr.att,1)
        corr_att(s,:) = squeeze(AAD_results_val_test_mix_with_ASR.corr.att(s,...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1) == params.tlag_min,...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s) == params.lambdas,:));
        corr_ign(s,:) = squeeze(AAD_results_val_test_mix_with_ASR.corr.ign(s,...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1) == params.tlag_min,...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s) == params.lambdas,:));
    end
    corr_att_ign = corr_att - corr_ign;
    
    %% ISC_same - ISC_other
    ISC_same_score = sum(ISC_same_ind(1:3,:),1);
    ISC_other_score = sum(ISC_other_ind(1:3,:),1);
    ISC_same_other = ISC_same_score - ISC_other_score;
    
    %% Mean SpE
    entropy_chan_seg_avg = mean(mean(SpE_results.entropy_norm,3),2);
    

%% Plotting
    figure('Units','centimeters','Position',[22,5,18,8]);

    %% AAD ~ SpE
    [rho_att_entr,p_att_entr] = corr(mean(corr_att_ign,2),mean(mean(SpE_results.entropy_norm,3),2),...
        'type','Spearman');
    
    sp1 = subplot(1,2,1);
    scatter(mean(corr_att_ign,2),mean(mean(SpE_results.entropy_norm,3),2),5,...
            'MarkerEdgeColor','none','MarkerFaceColor','k');
    xlabel('Attentional Gain (Corr_{att}-Corr_{ign}) [a.u.]');
    ylabel('Spectral Entropy [a.u.]');
    xlim([-0.005,0.06]);
    ylim([0.74,1]);
    lsline(sp1);
    
    %% ISC ~ SpE
    [rho_ISC_entr, p_ISC_entr] = corr(ISC_same_other',entropy_chan_seg_avg,...
        'type','Spearman');
    
    sp2 = subplot(1,2,2);
    scatter(ISC_same_other',entropy_chan_seg_avg,5,...
            'MarkerEdgeColor','none','MarkerFaceColor','k');
    xlabel('Attentional Effect (ISC_{same} - ISC_{other}) [a.u.]');
    ylabel('Spectral Entropy [a.u.]');
    xlim([-0.01,0.02]);
    ylim([0.74,1]);
    lsline(sp2);
    
    annotation('textbox',[0.06,0.88,0.1,0.1],'String','A','FontSize',10,'EdgeColor','none','FontWeight','bold');
    annotation('textbox',[0.5,0.88,0.1,0.1],'String','B','FontSize',10,'EdgeColor','none','FontWeight','bold');
    
    saveas(gcf, [fileout_plot_AAD_ISC_SpE_supp_figure,'.jpg']);
    saveas(gcf, [fileout_plot_AAD_ISC_SpE_supp_figure,'.png']);
    saveas(gcf, [fileout_plot_AAD_ISC_SpE_supp_figure,'.svg']);
    saveas(gcf, [fileout_plot_AAD_ISC_SpE_supp_figure,'.eps'],'epsc');
    saveas(gcf, [fileout_plot_AAD_ISC_SpE_supp_figure,'.tiff']);
    saveas(gcf, [fileout_plot_AAD_ISC_SpE_supp_figure,'.pdf']);
    close;



end

