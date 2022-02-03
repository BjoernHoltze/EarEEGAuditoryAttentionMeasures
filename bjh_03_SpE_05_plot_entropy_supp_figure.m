function bjh_03_SpE_05_plot_entropy_supp_figure(fileout_SpE_plot_supp_figure,SpE_results)
%% plots spectra and their entropy over time for different EEG channels
% 
% input:    SpE_results.freq_vec:       frequency vector of computed spectrum 
%           SpE_results.psd_norm:       normed power spectral density
%                                       (subj x seg x freq x chan)
%           SpE_results.entropy_norm:   normed spectral entropy 
%                                       (subj x seg x chan)
% author: Björn Holtze
% date: 14.01.22

    nb_seg = size(SpE_results.entropy_norm,2);

    figure('units','centimeters','outerposition',[0 0 18 25])
    trial_no = 1:nb_seg;
    for subj = 1:size(SpE_results.entropy_norm_chan_avg,1)
        % calculate statistics
        [rho_entr_time,p_entr_time] = corr(trial_no',squeeze(mean(SpE_results.entropy_norm(subj,:,:),3))',...
            'type','Spearman');
        % plot spectral entropy over time
        sp = subplot(6,6,subj);
        scatter(1:nb_seg,squeeze(mean(SpE_results.entropy_norm(subj,:,:),3)),10,'.');
        xlim([0,nb_seg+0.5]);
        ylim([0.65,1]);
        xline(nb_seg/3);
        xline(nb_seg*(2/3));
        if subj == 1
            xlabel('Time [min].');
            ylabel('Normed Entropy');
            sp.XLabel.FontSize = 7;
            sp.YLabel.FontSize = 7;
        end
        sp.XTick = [0,10,20,30];
        sp.XAxis.FontSize = 6;
        sp.YAxis.FontSize = 6;
        p_title = significant_stars(p_entr_time);
        t_txt = title(['rho: ',num2str(round(rho_entr_time,2)),' ',p_title]);
        t_txt.FontWeight = 'normal';
        t_txt.FontSize = 7;
    end
    saveas(gcf,[fileout_SpE_plot_supp_figure,'.jpg']);
    saveas(gcf,[fileout_SpE_plot_supp_figure,'.png']);
    saveas(gcf,[fileout_SpE_plot_supp_figure,'.svg']);
    saveas(gcf,[fileout_SpE_plot_supp_figure,'.eps'],'epsc');
    saveas(gcf,[fileout_SpE_plot_supp_figure,'.tiff']);
    saveas(gcf,[fileout_SpE_plot_supp_figure,'.pdf']);
    close;
    
   
end


