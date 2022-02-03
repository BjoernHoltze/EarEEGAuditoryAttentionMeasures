function bjh_01_AAD_07_plot_AAD_figure(file_out_plot_AAD_figure,AAD_results_val_test_mix_with_ASR,...
    AAD_results_val_test_mix_without_ASR,AAD_results_val_test_sep_with_ASR,params)

%% define figure parameters
figure('Units','centimeters','Position',[2,2,18,17]);
box_color = [0.8,0.8,0.8];

%% with vs. without ASR    
    % mutually-chosen hyperparameters
    % Wilcox signed rank test
    [p_ASR_ga,~,~] = signrank(AAD_results_val_test_mix_without_ASR.acc.ga,...
        AAD_results_val_test_mix_with_ASR.acc.ga,'method','approximate');

    [~,chance_level,~] = check_classifier(30/2,0.05,2);
    chance_level_per = chance_level*100;
    
    sp_1 = subplot(2,2,1);
    boxplot(cat(2,AAD_results_val_test_mix_without_ASR.acc.ga,AAD_results_val_test_mix_with_ASR.acc.ga),...
        'Colors',[0.5,0.5,0.5],'Symbol','.');
    hold on;
    set(findobj(sp_1,'type','line'),'LineWidth',1);
    boxes = findobj(sp_1,'Tag','Box');
    box_col =  cat(1,box_color,box_color);
    for b = 1:size(boxes,1)
        patch(get(boxes(b),'XData'),get(boxes(b),'YData'),box_col(b,:),...
            'EdgeColor','none','FaceAlpha',1);
    end
    line(cat(2,ones(size(AAD_results_val_test_mix_without_ASR.acc.ga,1),1),...
        ones(size(AAD_results_val_test_mix_with_ASR.acc.ga,1),1)+1)',...
        cat(2,AAD_results_val_test_mix_without_ASR.acc.ga,AAD_results_val_test_mix_with_ASR.acc.ga)',...
        'color','k','linestyle',':');
    line([-0.5,2.5],[chance_level_per,chance_level_per],'Color',[0.5,0.5,0.5]);
    sp_1.Children = sp_1.Children([end,end-2:end-1,1:end-3,1]);
    set(gca, 'Layer', 'top');
    sp_1.XTickLabel = {'without','with'};
    sp_1.XLabel.String = 'Artifact Correction';
    sp_1.XLabel.FontSize = 10;
    sp_1.YLabel.String = 'Decoding Accuracy [%]';
    sp_1.YLabel.FontSize =10;
    ylim([30,105]);
    s_ASR_ga = significant_stars(p_ASR_ga);
    sig_txt_ASR_ga = text(1.5,101,s_ASR_ga,'FontSize',9);
    sig_txt_ASR_ga.HorizontalAlignment = 'center';
    sp_1.Position(1) = 0.11;
    sp_1.Box = 'off';
    title('Effect of Artifact Correction');
    sp_1.Title.Position(2) = 107;
    
%% hyperparameter selection
    sp_2 = subplot(2,2,2);
    imagesc((params.tlag_min+params.tlag_max)/2,log10(params.lambdas),...
        squeeze(mean(AAD_results_val_test_mix_with_ASR.acc.cv,1))');
    caxis([35,100]);
    set(gca,'YDir','normal');
    cb = colorbar;
    cmap = colormap('jet');
    xlabel('Time Lag [ms]');
    ylabel('Log_{10} Reg. Parameter [a.u.]');
    sp_2.XLabel.FontSize = 10;
    sp_2.YLabel.FontSize = 10;
    cb.Label.String = 'Decoding Accuracy [%]';
    cb.Label.FontSize = 10;
    cb.Position(1) = 0.86;
    cb.Position([2,4]) = sp_1.Position([2,4]);
    cb.Position(3) = 0.015;
    hold on;
    % mutually-chosen hyperparameters
    rectangle('Position',[mean(AAD_results_val_test_mix_with_ASR.optimal_params.ga_tlag(1)+...
        AAD_results_val_test_mix_with_ASR.optimal_params.ga_tlag(2))/2-7.5,...
        log10(AAD_results_val_test_mix_with_ASR.optimal_params.ga_lambda)-0.3,15,0.6],...
        'EdgeColor','k','LineWidth',1.5);
    % individually-chosen hyperparameters
    m = length(cmap);
    % check whether individual optimal hyperparameters are the same for multiple subjects
    same_opt_params = cell(1,size(AAD_results_val_test_mix_with_ASR.acc.cv,1));
    for s = 1:size(AAD_results_val_test_mix_with_ASR.acc.cv,1)
        same_opt_params{s} = find(((AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1) == ...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(:,1)) & ...
            (AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s) ==...
            AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda)));
    end
    for s = 1:size(AAD_results_val_test_mix_with_ASR.acc.cv,1)
        if numel(same_opt_params{s}) == 1
            index = fix((AAD_results_val_test_mix_with_ASR.acc.ind(s)-cb.Limits(1))...
                /(cb.Limits(2)-cb.Limits(1))*m)+1;
            RGB = ind2rgb(index,cmap);
            plot((AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1)+AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,2))/2,...
                log10(AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s)),'o','MarkerSize',4,... % equation y = 0.1x - 2
                'MarkerEdgeColor','k','MarkerFaceColor',RGB);
        elseif numel(same_opt_params{s}) == 2
            if min(same_opt_params{s}) == s
                plot((AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1)+AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,2))/2,...
                    log10(AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s))-0.15,'o','MarkerSize',4,... % equation y = 0.1x - 2
                    'MarkerEdgeColor','k','MarkerFaceColor',RGB);
            elseif max(same_opt_params{s}) == s
                plot((AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1)+AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,2))/2,...
                    log10(AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s))+0.15,'o','MarkerSize',4,... % equation y = 0.1x - 2
                    'MarkerEdgeColor','k','MarkerFaceColor',RGB);
            end
        elseif numel(same_opt_params{s}) == 3
            if min(same_opt_params{s}) == s
                plot((AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1)+AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,2))/2,...
                    log10(AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s))-0.3,'o','MarkerSize',4,... % equation y = 0.1x - 2
                    'MarkerEdgeColor','k','MarkerFaceColor',RGB);
            elseif max(same_opt_params{s}) == s
                plot((AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1)+AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,2))/2,...
                    log10(AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s))+0.3,'o','MarkerSize',4,... % equation y = 0.1x - 2
                    'MarkerEdgeColor','k','MarkerFaceColor',RGB);
            else
                plot((AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,1)+AAD_results_val_test_mix_with_ASR.optimal_params.ind_tlag(s,2))/2,...
                    log10(AAD_results_val_test_mix_with_ASR.optimal_params.ind_lambda(s)),'o','MarkerSize',4,... % equation y = 0.1x - 2
                    'MarkerEdgeColor','k','MarkerFaceColor',RGB);
            end
        else
            error('More than 3 subjects have the same set of optimal parameteres');
        end
    end
    sp_2.Position(1) = 0.52;
    sp_2.Position(2:4) = sp_1.Position(2:4);
    sp_2.YTick = -5:1:5;
    sp_2.Box = 'off';
    title('Optimal Decoding Hyperparameters');
    sp_2.Title.Position(2) = 5.8;
    
%% individual vs. mutually-chosen hyperparameter (standard cross-validation)
    % Wilcox signed rank test
    [p_ind_ga_mix,~,~] = signrank(AAD_results_val_test_mix_with_ASR.acc.ga,...
        AAD_results_val_test_mix_with_ASR.acc.ind,'method','approximate');

    [~,chance_level,~] = check_classifier(30/2,0.05,2);
    chance_level_per = chance_level*100;
    
    sp_3 = subplot(2,2,3);
    boxplot(cat(2,AAD_results_val_test_mix_with_ASR.acc.ga,AAD_results_val_test_mix_with_ASR.acc.ind),'Colors',...
        [0.5,0.5,0.5],'Symbol','.');
    hold on;
    set(findobj(sp_3,'type','line'),'LineWidth',1);
    boxes = findobj(sp_3,'Tag','Box');
    box_col =  cat(1,box_color,box_color);
    for b = 1:size(boxes,1)
        patch(get(boxes(b),'XData'),get(boxes(b),'YData'),box_col(b,:),...
            'EdgeColor','none','FaceAlpha',1);
    end
    line(cat(2,ones(size(AAD_results_val_test_mix_with_ASR.acc.ga,1),1),...
        ones(size(AAD_results_val_test_mix_with_ASR.acc.ind,1),1)+1)',...
        cat(2,AAD_results_val_test_mix_with_ASR.acc.ga,AAD_results_val_test_mix_with_ASR.acc.ind)',...
        'color','k','linestyle',':');
    line([-0.5,2.5],[chance_level_per,chance_level_per],'Color',[0.5,0.5,0.5]);
    sp_3.Children = sp_3.Children([end,end-2:end-1,1:end-3,1]);
    set(gca, 'Layer', 'top');
    sp_3.XTickLabel = {'group-level','individual'};
    sp_3.XLabel.String = 'Hyperparameters';
    sp_3.XLabel.FontSize = 10;
    sp_3.YLabel.String = 'Decoding Accuracy [%]';
    sp_3.YLabel.FontSize = 10;
    ylim([30,105]);
    s_ind_ga_mix = significant_stars(p_ind_ga_mix);
    sig_txt_ind_ga_mix = text(1.5,101,s_ind_ga_mix,'FontSize',9);
    sig_txt_ind_ga_mix.HorizontalAlignment = 'center';
    sp_3.Position([1,3,4]) = sp_1.Position([1,3,4]);
    sp_3.Position(2) = 0.11;
    title(['Effect of Individual Hyperparameters', newline,'Standard Cross-Validation']);
    sp_3.Box = 'off';
    
%% individual vs. mutually-chosen hyperparameter (nested cross-validation)
    % Wilcox signed rank test
    [p_ind_ga_sep,~,~] = signrank(mean(AAD_results_val_test_sep_with_ASR.acc.ga,2),...
        mean(AAD_results_val_test_sep_with_ASR.acc.ind,2),'method','approximate');

    [~,chance_level,~] = check_classifier(10/2,0.05,2);
    chance_level_per = chance_level*100;
    
    sp_4 = subplot(2,2,4);
    boxplot(cat(2,mean(AAD_results_val_test_sep_with_ASR.acc.ga,2),...
        mean(AAD_results_val_test_sep_with_ASR.acc.ind,2)),'Colors',...
        [0.5,0.5,0.5],'Symbol','.');
    hold on;
    set(findobj(sp_4,'type','line'),'LineWidth',1);
    boxes = findobj(sp_4,'Tag','Box');
    box_col =  cat(1,box_color,box_color);
    for b = 1:size(boxes,1)
        patch(get(boxes(b),'XData'),get(boxes(b),'YData'),box_col(b,:),...
            'EdgeColor','none','FaceAlpha',1);
    end
    line(cat(2,ones(size(AAD_results_val_test_sep_with_ASR.acc.ga,1),1),...
        ones(size(AAD_results_val_test_sep_with_ASR.acc.ind,1),1)+1)',...
        cat(2,mean(AAD_results_val_test_sep_with_ASR.acc.ga,2),...
        mean(AAD_results_val_test_sep_with_ASR.acc.ind,2))',...
        'color','k','linestyle',':');
    line([-0.5,2.5],[chance_level_per,chance_level_per],'Color',[0.5,0.5,0.5]);
    sp_4.Children = sp_4.Children([end,end-2:end-1,1:end-3,1]);
    set(gca, 'Layer', 'top');
    sp_4.XTickLabel = {'group-level','individual'};
    sp_4.XLabel.String = 'Hyperparameters';
    sp_4.XLabel.FontSize = 10;
    sp_4.YLabel.String = 'Decoding Accuracy [%]';
    sp_4.YLabel.FontSize = 10;
    ylim([30,105]);
    s_ind_ga_sep = significant_stars(p_ind_ga_sep);
    sig_txt_ind_ga_sep = text(1.5,101,s_ind_ga_sep,'FontSize',9);
    sig_txt_ind_ga_sep.HorizontalAlignment = 'center';
    sp_4.Position(1) = sp_2.Position(1);
    sp_4.Position([2,3,4]) = sp_3.Position([2,3,4]);
    title(['Effect of Individual Hyperparameters', newline,'Nested Cross-Validation']);
    sp_4.Box = 'off';
    
    a = annotation('textbox',[0.044,0.86,0.1,0.1],'String','A','FontSize',10,'EdgeColor','none','FontWeight','bold');
    b = annotation('textbox',[0.46,0.86,0.1,0.1],'String','B','FontSize',10,'EdgeColor','none','FontWeight','bold');
    c = annotation('textbox',[0.044,0.403,0.1,0.1],'String','C','FontSize',10,'EdgeColor','none','FontWeight','bold');
    d = annotation('textbox',[0.46,0.403,0.1,0.1],'String','D','FontSize',10,'EdgeColor','none','FontWeight','bold');
    pause(0.5);
    a.Position(2) = 0.925;
    b.Position(2) = 0.925;
    c.Position(2) = 0.447;
    d.Position(2) = 0.447;
    
    saveas(gcf,[file_out_plot_AAD_figure,'.jpg']);
    saveas(gcf,[file_out_plot_AAD_figure,'.png']);
    saveas(gcf,[file_out_plot_AAD_figure,'.svg']);
    saveas(gcf,[file_out_plot_AAD_figure,'.eps'],'epsc');
    saveas(gcf,[file_out_plot_AAD_figure,'.tiff']);
    saveas(gcf,[file_out_plot_AAD_figure,'.pdf']);
    close;
    
end

