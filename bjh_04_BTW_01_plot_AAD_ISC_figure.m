function bjh_04_BTW_01_plot_AAD_ISC_figure(fileout_BTW_plot_AAD_ISC_figure,...
    ISC_same_ind,ISC_other_ind,AAD_results_val_test_mix_with_ASR,params)
    
    % compute ISC scores
    ISC_same_score = sum(ISC_same_ind(1:3,:),1);
    ISC_other_score = sum(ISC_other_ind(1:3,:),1);
    ISC_same_other = ISC_same_score - ISC_other_score;
    
    % compute corr_att-corr_ign
    corr_att = zeros(1,36);
    corr_ign = zeros(1,36);
    for s = 1:size(ISC_same_score,2)
        corr_att(s) = mean(AAD_results_val_test_mix_with_ASR.corr.att(s,...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1) == params.tlag_min,...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s) == params.lambdas,:),4);
        corr_ign(s) = mean(AAD_results_val_test_mix_with_ASR.corr.ign(s,...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1) == params.tlag_min,...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s) == params.lambdas,:),4);
    end
    corr_att_ign = corr_att - corr_ign; 
    
    % statistics
    [rho_AAD_ISC_diff,p_AAD_ISC_diff] = corr(corr_att_ign',ISC_same_other',....
        'type','Spearman','tail','right');

    % Plotting
    figure('units','centimeters','outerposition',[10 2 8.5 11])
    % backward model accuracy ~ ISC_same score
    AAD_ISC_diff = axes;
    scatter(corr_att_ign',ISC_same_other',5,...
        'MarkerEdgeColor','none','MarkerFaceColor','k');
    xlim([-0.005,0.06]);
    lsline(AAD_ISC_diff);
    AAD_ISC_diff.YAxis.Exponent = -3;
    AAD_ISC_diff.XAxis.Exponent = -2;
    xlabel('Attentional Gain (Corr_{att}-Corr_{ign})');
    ylabel('Attentional Effect (ISC_{same} - ISC_{other})');
    title(['Relation between Attentional Effects',newline,...
        'of Speech Envelope Tracking',newline,'and ISC']);
    AAD_ISC_diff.XAxis.Label.Position(1) = 0.025;

    saveas(gcf, [fileout_BTW_plot_AAD_ISC_figure,'.jpg']);
    saveas(gcf, [fileout_BTW_plot_AAD_ISC_figure,'.png']);
    saveas(gcf, [fileout_BTW_plot_AAD_ISC_figure,'.svg']);
    saveas(gcf, [fileout_BTW_plot_AAD_ISC_figure,'.eps'],'epsc');
    saveas(gcf, [fileout_BTW_plot_AAD_ISC_figure,'.tiff']);
    saveas(gcf, [fileout_BTW_plot_AAD_ISC_figure,'.pdf']);
    close;
end

