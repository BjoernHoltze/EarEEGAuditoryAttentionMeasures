function bjh_02_ISC_05_plot_ISC_figure(CONFIGPATH,fileout_plot_ISC_figure,ISC_same_ind,ISC_other_ind,ALLEEG_ISC_filt,...
                ISC_same_ind_chance,ISC_other_ind_chance,ISC_left_ind,ISC_right_ind,A_cond_indep,params)

% load subject information (id, attended channel, selected blocks)
load([CONFIGPATH,'subj_info.mat'],'subj_info');

attended_ch = [subj_info.attended_ch];

%% define function parameters
    figure('Units','centimeters','Position',[5,2,18,19]);

%% ISc same/other
    ISC_same_score = sum(ISC_same_ind(1:params.Ncomp,:),1);
    ISC_other_score = sum(ISC_other_ind(1:params.Ncomp,:),1);

    pure_red_0_5 = [236,168,139]/255;
    pure_blue_0_5 = [127,184,222]/255; 
    pure_red = [0.8500, 0.3250, 0.0980];  
    pure_blue = [0,0.4470,0.7410]; 
    
    % statistics
        % check for normality
        kstest(ISC_same_score-ISC_other_score);
        % paired-sample t-test
        [~,p_same_other] = ttest(ISC_same_score,ISC_other_score);

    % plot boxplots
    sp_1 = subplot(4,6,[1,2,3,7,8,9]); 
    boxplot([ISC_same_score;ISC_other_score]','Colors',...
        [pure_red;pure_blue]);  
    set(findobj(sp_1,'type','line'),'LineWidth',1);
    boxes = findobj(sp_1,'Tag','Box');
    box_col = [pure_blue_0_5;pure_red_0_5];
    for b = 1:size(boxes,1)
    patch(get(boxes(b),'XData'),get(boxes(b),'YData'),box_col(b,:),'EdgeColor','none','FaceAlpha',1);
    end
    hold on;
    % connect individual data points by line
    line([ones(size(ISC_same_score,2),1),ones(size(ISC_same_score,2),1)+1]',...
        [ISC_same_score',ISC_other_score']','color','k','linestyle',':');
    % draw in chance level as a line
    ISC_same_chance_scores = squeeze(sum(ISC_same_ind_chance(1:params.Ncomp,:,:),1));
    ISC_same_chance_scores_95 = prctile(ISC_same_chance_scores(:),95);
    line([0.5,1.5],[ISC_same_chance_scores_95,ISC_same_chance_scores_95],'Color',[0.5,0.5,0.5]);
    ISC_other_chance_scores = squeeze(sum(ISC_other_ind_chance(1:3,:,:),1));
    ISC_other_chance_scores_95 = prctile(ISC_other_chance_scores(:),95);
    line([1.5,2.5],[ISC_other_chance_scores_95,ISC_other_chance_scores_95],'Color',[0.5,0.5,0.5]);
    set(sp_1, 'Children',sp_1.Children([end,end-1,end-2,1:end-3]));
    sp_1.XTickLabel = {'Same','Other'};
    sp_1.XLabel.String = 'Attended Story';
    sp_1.XLabel.FontSize = 10;
    sp_1.YLabel.String = 'ISC Sum Score [a.u.]';
    sp_1.YLabel.FontSize = 10;
    sp_1.YAxis.Exponent = -3;
    h = findobj(gcf,'tag','Outliers');
    set(h,'Marker','o','MarkerSize',2,'MarkerEdgeColor',pure_red,'MarkerFaceColor',pure_red);
    star_same_other = significant_stars(p_same_other);
    star_same_other_txt = text(1.5,0.015,star_same_other);
    star_same_other_txt.HorizontalAlignment = 'center';
    sp_1.Position(1) = 0.09;
    sp_1.Position(2) = 0.57;
    sp_1.Box = 'off';
    sp_1.Title.String = ['Attentional Effect on ISC',newline,...
        'Same vs. Other'];
    sp_1.Title.FontSize = 10;
   
    
%% Scatter Plot (ISC_left ~ ISC_right)
    ISC_left_ind_score = sum(ISC_left_ind(1:params.Ncomp,:),1);
    ISC_right_ind_score = sum(ISC_right_ind(1:params.Ncomp,:),1);
    
    pure_yell = [0.9290 0.6940 0.1250];
    pure_lila = [0.4940 0.1840 0.5560];
    
    sp_2 = subplot(4,6,[4,5,6,10,11,12]);
    plot(ISC_left_ind_score(attended_ch==1),ISC_right_ind_score(attended_ch==1),'o',...
        'MarkerFaceColor',pure_lila,'MarkerEdgeColor',pure_lila,'MarkerSize',4);
    hold on; 
    plot(ISC_left_ind_score(attended_ch==2),ISC_right_ind_score(attended_ch==2),'o',...
        'MarkerFaceColor',pure_yell,'MarkerEdgeColor',pure_yell,'MarkerSize',4);
    hold off;
    sp_2.XAxis.Limits = [-0.006,0.02];
    sp_2.XAxis.Exponent = -3;
    sp_2.YAxis.Exponent = -3;
    legend('attended left','attended right');
    sp_2.Legend.FontSize = 8;
    sp_2.Legend.Position = [0.769,0.85,0.1795,0.0441];
    
    xlabel('ISC_{left} Sum Score [a.u.]');
    ylabel('ISC_{right} Sum Score [a.u.]');
    sp_2.Position(2:4) = sp_1.Position(2:4);
    sp_2.Position(1) = 0.575;
    sp_2.Box = 'off';
    sp_2.XLabel.FontSize = 10;
    sp_2.YLabel.FontSize = 10;
    sp_2.Title.String = ['Attentional Effect on ISC',newline,...
        'Left vs. Right'];
    sp_2.Title.FontSize = 10;


%% Plot ISC_same/ISC_other - Individual Components
    
    ISC_comp = zeros(params.Ncomp,2);
    ISC_comp_chance = zeros(params.Ncomp,2);
    h = zeros(1,3);
    p_comp = zeros(1,3);
    for c = 1:params.Ncomp
        % create matrix for bar plot
        ISC_comp(c,1) = mean(ISC_same_ind(c,:));
        ISC_comp(c,2) = mean(ISC_other_ind(c,:));
        ISC_comp_chance(c,1) = prctile(ISC_same_ind_chance(c,:),95);
        ISC_comp_chance(c,2) = prctile(ISC_other_ind_chance(c,:),95);  
        % check for normality
        h(c) = kstest(ISC_same_ind(c,:)-ISC_other_ind(c,:));
        % paired-sample t-test
        [~,p_comp(c)] = ttest(ISC_same_ind(c,:),ISC_other_ind(c,:));
    end        
    
    % plot bar plot
    sp_3 = subplot(4,6,13:18);
    b = bar(ISC_comp,0.4);
    b(1).EdgeColor = pure_red;
    b(1).FaceColor = pure_red_0_5;
    b(2).EdgeColor = pure_blue;
    b(2).FaceColor = pure_blue_0_5;
    sp_3.YLabel.String = 'ISC Score [a.u.]';   
    sp_3.YLabel.FontSize = 10;
    sp_3.XTickLabel = {};
    sp_3.XAxis.TickLabels = {'Component 1','Component 2','Component 3'};
    sp_3.XAxis.FontSize = 10;
    
    % plot chance level
    hold on; 
    bar(ISC_comp_chance,0.4,'FaceColor','k','EdgeColor','k',...
        'FaceAlpha',0.2,'EdgeAlpha',0.2);
    hold off;
    legend('attended same','attended other');
    sp_3.Legend.FontSize = 8;
    sp_3.Legend.Position = [0.769,0.416,0.1795,0.0441];
    
    ylim([0,0.007]);
    
    % plot significance
    for c = 1:params.Ncomp
        star_comp = significant_stars(p_comp(c));
        if c == 3
            star_comp_txt = text(c,0.0035,star_comp);
        else
            star_comp_txt = text(c,0.005,star_comp);
        end
        star_comp_txt.HorizontalAlignment = 'center';
        star_comp_txt.FontSize = 8;
    end
    sp_3.Box = 'off';
    sp_3.Position(1) = sp_1.Position(1);
    sp_3.Position(3) = (sp_2.Position(1)+sp_2.Position(3))-sp_1.Position(1);
    sp_3.Position(2) = 0.312;
    sp_3.Position(4) = 0.13;
    sp_3.Title.String = 'ISC Scores of Individual Components';
    sp_3.Title.FontSize = 10;
    sp_3.Title.Position(2) = 0.007;
    
    
    
%% plot forward model
    isc_fmodel.setname = '';
    isc_fmodel.times = 1;
    cEEGrid_chanlabels = {ALLEEG_ISC_filt(1).chanlocs.labels};
    
    % Component 1
    sp_4 = subplot(4,6,[19,20]);
    isc_fmodel.data = A_cond_indep(:,1);
    cfg.colormap = 'jet';
    cEEGrid_topoplot(isc_fmodel,cEEGrid_chanlabels,cfg,1,'interpolation','v4',...
        'newfig','no','comment','no','zlim',[-0.3,0.3]);
    cb = colorbar;
    cb.Visible = 'off';
    pause(0.1);
    sp_4.Position(1) = 0.14;
    sp_4.Position(2) = 0.04;
    for child = 1:size(sp_4.Children,1)
        sp_4.Children(child).LineWidth = 1;
    end
    t_4 = title('Component 1');
    t_4.Position(2) = -110;
    t_4.FontWeight = 'Normal';
    t_4.FontSize = 10;
    
    % Component 2
    sp_5 = subplot(4,6,[21,22]);
    isc_fmodel.data = A_cond_indep(:,2);
    cfg.colormap = 'jet';
    cEEGrid_topoplot(isc_fmodel,cEEGrid_chanlabels,cfg,1,'interpolation','v4',...
        'newfig','no','comment','no','zlim',[-0.3,0.3]);
    cb = colorbar;
    cb.Visible = 'off';
    pause(0.1);
    sp_5.Position(1) = 0.435;
    sp_5.Position(2) = 0.04;
    for child = 1:size(sp_5.Children,1)
        sp_5.Children(child).LineWidth = 1;
    end
    t_5 = title('Component 2');
    t_5.Position(2) = -110;
    t_5.FontWeight = 'Normal';
    t_5.FontSize = 10;    
    
    % Component 3
    sp_6 = subplot(4,6,[23,24]);
    isc_fmodel.data = A_cond_indep(:,3);
    cfg.colormap = 'jet';
    cEEGrid_topoplot(isc_fmodel,cEEGrid_chanlabels,cfg,1,'interpolation','v4',...
        'newfig','no','comment','no','zlim',[-0.3,0.3]);
    cb = colorbar;
    pause(0.1);
    sp_6.Position(1) = 0.732;
    sp_6.Position(2) = 0.04;
    for child = 1:size(sp_6.Children,1)
        sp_6.Children(child).LineWidth = 1;
    end
    cb.Position(1) = 0.905;
    cb.Label.String = '[a.u.]';
    cb.Label.FontSize = 8;
    t_6 = title('Component 3');
    t_6.Position(2) = -110;
    t_6.FontWeight = 'Normal';
    t_6.FontSize = 10; 
    sp_title = annotation('textbox',[0.45,0.19,0.1,0.1],'String','Spatial Patterns of Individual Components',...
        'FontSize',10,'FontWeight','bold','EdgeColor','none',...
        'HorizontalAlignment','center'); 

    a = annotation('textbox',[0.035,0.933,0.1,0.1],'String','A','FontSize',10,'EdgeColor','none','FontWeight','bold');
    b = annotation('textbox',[0.5,0.933,0.1,0.1],'String','B','FontSize',10,'EdgeColor','none','FontWeight','bold');
    c = annotation('textbox',[0.035,0.4590,0.1,0.1],'String','C','FontSize',10,'EdgeColor','none','FontWeight','bold');
    d = annotation('textbox',[0.035,0.19,0.1,0.1],'String','D','FontSize',10,'EdgeColor','none','FontWeight','bold');
    pause(0.5);
    a.Position(2) = 0.937;
    b.Position(2) = 0.937;
    c.Position(2) = 0.46;
    d.Position(2) = 0.195;
    sp_title.Position(1) = 0.308;
    sp_title.Position(2) = 0.195;
    
    saveas(gcf,[fileout_plot_ISC_figure,'.jpg']);
    saveas(gcf,[fileout_plot_ISC_figure,'.png']);
    saveas(gcf,[fileout_plot_ISC_figure,'.svg']);
    saveas(gcf,[fileout_plot_ISC_figure,'.eps'],'epsc');
    saveas(gcf,[fileout_plot_ISC_figure,'.tiff']);
    saveas(gcf,[fileout_plot_ISC_figure,'.pdf']);
    close;

end

