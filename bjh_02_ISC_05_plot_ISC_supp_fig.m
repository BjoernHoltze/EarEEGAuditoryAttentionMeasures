function bjh_02_ISC_05_plot_ISC_supp_fig(file_out_ISC_forward_model_with_asr,ALLEEG_ISC_filt,...
    A_left,A_right,params)
%% plots forward model topographies
% 
% input:    file_out_plot_cond_indep_ISC: full_path including file name but 
%                                       without extension where
%                                       figure should bes stored.
%           A_cond_indep:       forward model projections
%           params:             [struct] structure containing number of
%                               components to be plotted
% 
% authore: Björn Holtze 
% date: 06.08.21

%% forward model topographies and individual condition independent ISC values
    isc_fmodel.setname = '';
    isc_fmodel.times = 1;
    cEEGrid_chanlabels = {ALLEEG_ISC_filt(1).chanlocs.labels};

    figure('Units','centimeters','Position',[10,6,18,12]);
    for condition = 1:2
        for c = 1:params.Ncomp
            % Condition inepdentent forward model
            if condition == 1
                isc_fmodel.data = A_left(:,c);
            elseif condition == 2
                isc_fmodel.data = A_right(:,c);
            end
            cfg.colormap = 'jet';
            sp_pos = (condition*3-3)+c;
            sp(sp_pos) = subplot(2,params.Ncomp,(condition*3-3)+c);
            cEEGrid_topoplot(isc_fmodel,cEEGrid_chanlabels,cfg,1,'interpolation','v4',...
                'newfig','no','comment','no','zlim',[-0.3,0.3]);
            cb(sp_pos) = colorbar;
            if sp_pos ~= 6
            cb(sp_pos).Visible = 'off';
            end
            title(['Component ',num2str(c)]);
        end
    end
    
    for sp_pos = 1:6
        sp(sp_pos).Position(1) = sp(sp_pos).Position(1) - 0.03;
        if sp_pos <= 3
            sp(sp_pos).Position(2) = 0.4;
        elseif sp_pos > 3
            sp(sp_pos).Position(2) = -0.05;
        end
        sp(sp_pos).Position(3) = 0.22;
        sp(sp_pos).Position(4) = 0.5;
        sp(sp_pos).Title.Position(2) = 400;
    end
    
    txt_left = text(-30,1230,'Left Condition');
    txt_left.HorizontalAlignment = 'center';
    txt_left.FontWeight = 'bold';
    txt_right = text(-30,530,'Right Condition');
    txt_right.HorizontalAlignment = 'center';
    txt_right.FontWeight = 'bold';
    
    cb(6).Position = [0.9, 0.091, 0.015, 0.21];
    cb(6).Label.String = '[a.u.]';
    
    saveas(gcf,[file_out_ISC_forward_model_with_asr,'.jpg']);
    saveas(gcf,[file_out_ISC_forward_model_with_asr,'.png']);
    saveas(gcf,[file_out_ISC_forward_model_with_asr,'.svg']);
    saveas(gcf,[file_out_ISC_forward_model_with_asr,'.eps'],'epsc');
    saveas(gcf,[file_out_ISC_forward_model_with_asr,'.tiff']);
    saveas(gcf,[file_out_ISC_forward_model_with_asr,'.pdf']);

    close;

end

