function bjh_03_SpE_05_plot_entropy_figure(fileout_SpE_plot_figure,SpE_results,params)
%% plots spectra and their entropy over time for different EEG channels
% 
% input:    SpE_results.freq_vec:       frequency vector of computed spectrum 
%           SpE_results.psd_norm:       normed power spectral density
%                                       (subj x seg x freq x chan)
%           SpE_results.entropy_norm:   normed spectral entropy 
%                                       (subj x seg x chan)
% author: Björn Holtze
% date: 14.01.22

%% define figure parameters
    figure('units','centimeters','outerposition',[10 2 8.5 20])
    pure_blue = [0,0.4470,0.7410]; 

%% grand average time-frequency plot
    nb_seg = size(SpE_results.entropy_norm,2);
    trial_no = 1:nb_seg;
    % plot spectrum
    sp1 = subplot(3,1,1);
    xline(nb_seg/3);
    hold on;
    xline(nb_seg*(2/3));
    imagesc(1:nb_seg,params.entr_freq_min:params.entr_freq_max,...
        100*squeeze(mean(SpE_results.psd_norm_chan_avg(:,:,SpE_results.freq_vec >= params.entr_freq_min & ...
        SpE_results.freq_vec <= params.entr_freq_max),1))',[1,5.3]);
    cb = colorbar;
    cb.Position(3) = 0.02;
    cb.Label.String = 'Power [\muV^2]';
    cb.Label.FontSize = 9;
    cb.Label.Position(1) = 2.8;
    sp1.Position(1) = 0.15;
    sp1.Position(2) = 0.72;
    sp1.Position(3) = 0.7;
    cb.Position(2) = sp1.Position(2);
    cb.Position(1) = 0.86;
    set(gca,'YDir','normal');
    sp1.XAxis.Limits = [0.5,nb_seg+0.5];
    sp1.XAxis.Label.String = 'Time [min]';
    sp1.YAxis.Label.String = 'Frequency [Hz]';
    text(33,34,'x10^{-2}','FontSize',7.5);
    sp1.Box = 'off';
    sp1.Title.String = 'Grand Average cEEGrid Spectogram';
    pause(0.1);
   
    
%% grand average entropy over time
    gr_avg_entropy = squeeze(mean(mean(SpE_results.entropy_norm,3),1));
    [rho_gr_avg_entr_time,p_gr_avg_entr_time] = corr(trial_no',gr_avg_entropy',...
        'type','Spearman');
    sp2 = subplot(3,1,2);
    errorbar(trial_no,gr_avg_entropy,std(mean(SpE_results.entropy_norm,3),0,1)./...
        sqrt(size(SpE_results.subj_info,2)),'o','MarkerSize',3,...
        'MarkerFaceColor',pure_blue);
    sp2.XAxis.Label.String = 'Time [min]';
    sp2.YAxis.Label.String = 'Spectral Entropy [a.u.]';
    sp2.XAxis.Limits = [0.5,30.5];
    xline(nb_seg/3);
    xline(nb_seg*(2/3));
    sp2.Position(1) = sp1.Position(1);
    sp2.Position(2) = 0.397;
    sp2.Position(3) = sp1.Position(3);
    sp2.Box = 'off';
    sp1.YAxis.Label.Position(1) = sp2.YAxis.Label.Position(1);
    sp2.Title.String = 'Grand Average Spectral Entropy';
    
%% alpha power over time
    gr_avg_alpha = mean(mean(SpE_results.psd_norm_chan_avg(:,:,SpE_results.freq_vec >= 8 & ...
        SpE_results.freq_vec <= 12),3),1);
    [rho_gr_avg_alpha_time,p_gr_avg_alpha_time] = corr(trial_no',gr_avg_alpha',...
        'type','Spearman');
    sp3 = subplot(3,1,3);
    errorbar(trial_no,gr_avg_alpha,std(mean(SpE_results.psd_norm_chan_avg(:,:,...
        SpE_results.freq_vec >= 8 & SpE_results.freq_vec <= 12),3),0,1)./...
        sqrt(size(SpE_results.subj_info,2)),'o','MarkerSize',3,...
        'MarkerFaceColor',pure_blue);
    sp3.XAxis.Label.String = 'Time [min]';
    sp3.YAxis.Label.String = 'Alpha Power [\muV^2]';
    sp3.XAxis.Limits = [0.5,30.5];
    xline(nb_seg/3);
    xline(nb_seg*(2/3));
    sp3.Position(1) = sp1.Position(1);
    sp3.Position(2) = 0.07;
    sp3.Position(3) = sp1.Position(3);
    sp3.Box = 'off';
    sp3.YAxis.Label.Position(1) = sp2.YAxis.Label.Position(1);
    sp3.Title.String = 'Grand Average Alpha Power';
    sp3.YAxis.Exponent = -2;

    annotation('textbox',[0.015,0.88,0.1,0.1],'String','A','FontSize',10,'EdgeColor','none','FontWeight','bold');
    annotation('textbox',[0.015,0.565,0.1,0.1],'String','B','FontSize',10,'EdgeColor','none','FontWeight','bold');
    annotation('textbox',[0.015,0.245,0.1,0.1],'String','C','FontSize',10,'EdgeColor','none','FontWeight','bold');  
    
    saveas(gcf,[fileout_SpE_plot_figure,'.jpg']);
    saveas(gcf,[fileout_SpE_plot_figure,'.png']);
    saveas(gcf,[fileout_SpE_plot_figure,'.svg']);
    saveas(gcf,[fileout_SpE_plot_figure,'.eps'],'epsc');
    saveas(gcf,[fileout_SpE_plot_figure,'.tiff']);
    saveas(gcf,[fileout_SpE_plot_figure,'.pdf']);
    close;

end

