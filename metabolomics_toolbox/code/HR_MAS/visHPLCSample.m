function h = visHPLCSample(sample,varargin)

% MTJ 2021

% Defaults

% sample = ;
% cmap = ;
% ylabelL = ;
% ylabelR = ;
% titleStr = ;

% if ~exist('dataField','var')
%     dataField = 'data';
% end

if ~isempty(varargin)
    
    if any(contains(varargin,'save'))
        saveFig = true;
    end
    if any(contains(varargin,'plotSum'))
        plotSum = true;
    end    
end

    for s = 1:length(sample)
        
        % Plot the absorbance   
            figure, hold on
            imagesc(sample(s).timepoints,sample(s).wavelengths, sample(s).dataScaled')
               
                colormap('jet') % Convert this to wavelengths? 
                    
                set(gca,'ydir','normal')
                xlabel('Time (min)')
                    ylabel('Wavelength (nm)')
                    set(gca,'ylim',[min(sample(s).wavelengths),max(sample(s).wavelengths)])
%                 c = colorbar('Location','westoutside');
                c = colorbar('Location','eastoutside');
                    ylabel(c, 'Absorbance (log-scaled, > 0)')   
        
        if exist('plotSum','var')
            % Plot the summed absorbance on top in bold white
                yyaxis right
%                     plot(sample(s).timepoints,sum(sample(s).data,2),'Color','w','LineWidth',2) %,'DisplayName',['Sum of Intensities for ',sample(s).filename],'Interpreter','none')
%                     ylabel('Total Absorbance (raw)')

                    plot(sample(s).timepoints,sum(sample(s).dataScaled,2),'Color','w','LineWidth',2) %,'DisplayName',['Sum of Intensities for ',sample(s).filename],'Interpreter','none')
                    ylabel('Total Absorbance (scaled)')
        end
        
        % Title for both
            title(sample(s).name,'Interpreter','none')
            set(gca,'xlim',[min(sample(s).timepoints),max(sample(s).timepoints)])

        % Saving 
            if exist('saveFig','var')
               fig = gca;saveas(fig,[fig.Title.String],'pdf')
            end
    end
end