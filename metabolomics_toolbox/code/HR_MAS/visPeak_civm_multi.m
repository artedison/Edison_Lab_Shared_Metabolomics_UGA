        % Add Public toolbox
            addpath(genpath('/Users/mjudge/Edison_Lab_Shared_Metabolomics_UGA')) 
        % Remove Private toolbox
            rmpath(genpath('/Users/mjudge/Edison_lab_UGA'))
            pause(1),clc

%% Get the datasets

    cd('/Users/mjudge/Dropbox (Edison_Lab@UGA)/Projects/clock/CIVM_paper_2/Acetate_QAX_project/NMRdata/processed')

    experiments = {'civm_ncrassa_qax_02';
                'civm_ncrassa_qax_03';
                'CIVM_ncrassa_qax_06_pgm2_gluc_1';
                'CIVM_ncrassa_qax_08';
                'CIVM_ncrassa_qax_10';
                'CIVM_ncrassa_qax_11';
                'paper_Sample_4';
                'CIVM_ncrassa_qax_12'};
            
% Open the above files
%     for e = 1:length(experiments)
%         cd(experiments{e})
%         edit('makeDirs*.m')
%         cd ..
%     end
    
    datasets = struct();    
    
    for d = 1:length(experiments)
        datasets(d).name = experiments{d};
        
        cd(experiments{d})
            temp = dir('*short.mat');
            datasets(d).file = temp(1);
            
        datasets(d).dataStruct = open(datasets(d).file.name);
        
        datasets(d).structName = fields(datasets(d).dataStruct);                      % 
        datasets(d).dataStruct = datasets(d).dataStruct.(datasets(d).structName{:});  % re- assign the unpacked data
        
        cd ..
        
    end

    clear('temp','experiments','d','tempstruct','e')
    cd ..
    
%% Plot peak from each sample

        noise = 1.298578433674705e+10; % got this as output from running stackSpectra on the data
        timeLimit = 12;
        
   plotSpec(1).region = [5,5.6]; % G1P
        plotSpec(1).peakName = 'glucose-1-phoshphate';
        
   plotSpec(2).region = [1.1498    1.22]; % Ethanol
        plotSpec(2).peakName = 'ethanol';
        
   plotSpec(3).region = [3.0    3.3]; % 3.2 singlet
        plotSpec(3).peakName = 'Singlet at 3.2';
        
   plotSpec(4).region = [2.3    2.6]; % succinate
        plotSpec(4).peakName = 'succinate';

   plotSpec(5).region = [3    4.5]; % succinate
        plotSpec(5).peakName = 'myo-inositol region';

   plotSpec(6).region = [-0.5 10];
        plotSpec(6).peakName = 'Whole Spectrum';
        
   plotSpec(7).region = [-0.0496    1.2207];
        plotSpec(7).peakName = 'EtOH cp DSS';   
        
   plotSpec(8).region = [1.8353    1.8621];
        plotSpec(8).peakName = 'quinic acid (quant)'; 

   plotSpec(9).region = [5.211    5.23968];
        plotSpec(9).peakName = 'glucose a (quant)'; 

   plotSpec(10).region = [1.8353    1.8621];
        plotSpec(10).peakName = 'glucose b (quant)'; 
        
%% 
p = 2;
clear ax fig 
    maxTime = 40;
    plotSpec(p).plotRes = 30;
    plotSpec(p).horzshift = .002;
    plotSpec(p).vertshift = 0.8 * noise;   
    quinicData = [1,2,5];
    for d = 1:length(datasets)

%     for d = 8%quinicData

        timepoints = [datasets(d).dataStruct.smoothedData.timepoints];    
        maxInd = max(find(timepoints<maxTime));
        ppm = datasets(d).dataStruct.ppm;
        reginds = fillRegion(plotSpec(p).region,ppm);
        ppm = ppm(reginds);
        matrix = vertcat(datasets(d).dataStruct.smoothedData.data);
        matrix = matrix(:,reginds);
            
             % Make a Stack Plot of the spectra:

                    [datasets(d).dataStruct.plotInds,datasets(d).dataStruct.plotIndsCat] = calc_stackPlotInds({datasets(d).dataStruct.smoothedData.data},plotSpec(p).plotRes,maxInd);
                    
                [~,plotSpec(p).params] = stackSpectra(matrix,ppm,...
                            plotSpec(p).horzshift,... 
                            plotSpec(p).vertshift,...
                            [datasets(d).dataStruct.plotTitle;' (',...
                                    num2str(timepoints(max(datasets(d).dataStruct.plotIndsCat))),...
                            ' h)'],...
                             'colors',datasets(d).dataStruct.colorsSmoothed,...
                             'plotSubset',datasets(d).dataStruct.plotIndsCat,...
                             'noWhiteShapes',...
                             'timeVect',timepoints);
%                          legend 'off'

                 
%                  highlightROIs([1.8353    1.8621]',max(matrix(:)),'horzshift', plotSpec(p).params.horzshift,'vertshift',plotSpec(p).params.vertshift, 'numberOfSpectra', size(matrix,1),'extension',0.1)
                 qaquant_max{d} = max(matrix,[],2);
                 qaquant_min{d} = min(matrix,[],2);
                 
        ax(d) = gca;
        fig(d) = gcf;

%         set(gca,'xlim',plotSpec.xlims)
%         set(gca,'ylim',plotSpec.ylims)
        
    end
    
 linkaxes(ax);
 
 
%% Saving
cd('/Users/mjudge/Dropbox (Edison_Lab@UGA)/Projects/clock/CIVM_paper_2/Acetate_QAX_project/NMRdata/processed')
mkdir('multisample_plotting_byPeak')
cd('multisample_plotting_byPeak')
regionName = [num2str(plotSpec(p).region(1)),'-',num2str(plotSpec(p).region(2)),' ppm'];
mkdir(regionName);
cd(regionName)

plotSpec(p).xlims = get(gca,'xlim');
plotSpec(p).ylims = get(gca,'ylim');


    for i = 1:length(ax)   
        printCleanPDF(fig(i),[datasets(i).structName{:},...
                                ' - Region ',...
                                    num2str(plotSpec(p).region(1)),'-',num2str(plotSpec(p).region(2)),...
                                ' ppm']);
    end

    close all
    
%%
    for i = 1:length(sn)  
        for j = 1:length(bs)
            n =  i + (j-1) * length(sn); 
            % Make the ax optOB_out.resultsject
                ax(j,i) = subplot(length(bs),length(sn),n);
                    hold on
                    
                    %% Make a Stack Plot of the spectra:
        
                        matrix = vertcat(pgm2_starve_gluc.Xsmoothed.data);
                        currentppm = pgm2_starve_gluc.ppmcat;
                        pgm2_starve_gluc.plotRes = 50;

                            [pgm2_starve_gluc.plotInds,pgm2_starve_gluc.plotIndsCat] = calc_stackPlotInds({pgm2_starve_gluc.Xsmoothed.data},pgm2_starve_gluc.plotRes);

                        pgm2_starve_gluc.horzshift = .001;
                        pgm2_starve_gluc.vertshift = 0.3;

                        plotTitle = {'pgm2 starved,'; 'glucose pulsed,';' then starved'};        

                        stackSpectra(matrix,currentppm,pgm2_starve_gluc.horzshift,pgm2_starve_gluc.vertshift,plotTitle,...
                                     'colors',pgm2_starve_gluc.colorsSmoothed,...
                                     'autoVert',...
                                     'plotSubset',pgm2_starve_gluc.plotIndsCat)
            
        end
    end
    linkaxes(ax(:),'xy');
        
    % Add a title
        sttl = suptitle('Results for opt_bucket Parameter Optimization','interpreter','none');
    
    % Add one x label across the bottom
        suplabel('Chemical Shift (ppm)','x');
    

%% QA and total spectral intensity Quant
 figure,hold on
    for d = 1
        m = qaquant_max{d};
            m = m/max(m(:));
        plot([datasets(d).dataStruct.smoothedData.timepoints],m)
    end 
    legend({'WT','qa-X','qa-4'})
    title('QA Consumption in Mutants')
    xlabel('Time (h)')
    ylabel('Scaled QA Peak Height (Peak Max)')
%     ylabel('Total Spectral Intensity')
    set(gca,'FontSize',20)
 figure,hold on
    for d = [1,2,5]
        m = qaquant_max{d}-qaquant_min{d};
            m = m/max(m(:));
        plot([datasets(d).dataStruct.smoothedData.timepoints],m)
    end 
    legend({'WT','qa-X','qa-4'})
    title('QA Consumption in Mutants')
    xlabel('Time (h)')
    ylabel('Scaled QA Peak Height (Max - min)')
%     ylabel('Total Spectral Intensity')
    set(gca,'FontSize',20)

%% EtOH and total spectral intensity Quant
 figure,hold on
    for d = 8
        m = qaquant_max{d};
            m = m/max(m(:));
        plot([datasets(d).dataStruct.smoothedData.timepoints],m)
    end 
    legend({'pgm-2'})
    title('EtOH Long Starve pgm-2')
    xlabel('Time (h)')
    ylabel('Scaled EtOH Peak Height (Peak Max)')
%     ylabel('Total Spectral Intensity')
    set(gca,'FontSize',20)
 figure,hold on
    for d = 8
        m = qaquant_max{d}-qaquant_min{d};
            m = m/max(m(:));
        plot([datasets(d).dataStruct.smoothedData.timepoints],m)
    end 
    legend({'pgm-2'})
    title('EtOH Long Starve pgm-2')
    xlabel('Time (h)')
    ylabel('Scaled EtOH Peak Height (Max - min)')
%     ylabel('Total Spectral Intensity')
    set(gca,'FontSize',20)
    
    
%% QA and total spectral intensity Quant
    for d = 1:length(datasets)
        figure
        m = qaquant_max{d};
            m = m/max(m(:));
        plot([datasets(d).dataStruct.smoothedData.timepoints],m)
        title(['EtOH in ',datasets(d).structName],'Interpreter','none')
        xlabel('time (h)')
    end 
    