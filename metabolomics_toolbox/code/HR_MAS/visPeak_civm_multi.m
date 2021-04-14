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
                'paper_Sample_4'};
            
% Open the above files
    for e = 1:length(experiments)
        cd(experiments{e})
        edit('makeDirs*.m')
        cd ..
    end
    
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

    clear('temp','experiments','d','tempstruct')
    cd ..
    
%% Plot peak from each sample

        noise = 1.298578433674705e+10; % got this as output from running stackSpectra on the data
        timeLimit = 12;
        
   plotSpec(1).region = [5,5.6]; % G1P
        plotSpec(1).peakName = 'glucose-1-phoshphate';
        
   plotSpec(2).region = [1.1498    1.2630]; % Ethanol
        plotSpec(2).peakName = 'ethanol';
        
   plotSpec(3).region = [3.0    3.3]; % 3.2 singlet
        plotSpec(3).peakName = 'Singlet at 3.2';
        
   plotSpec(4).region = [2.3    2.6]; % succinate
        plotSpec(4).peakName = 'succinate';

%% 
p = 4;
clear ax fig 

%     plotSpec(p).plotRes = 50;
%     plotSpec(p).horzshift = .001;
%     plotSpec(p).vertshift = 0.2 * noise;   

    for d = 1:length(datasets)

        timepoints = [datasets(d).dataStruct.smoothedData.timepoints];        
        ppm = datasets(d).dataStruct.ppm;
        reginds = fillRegion(plotSpec(p).region,ppm);
        ppm = ppm;
        matrix = vertcat(datasets(d).dataStruct.smoothedData.data);
        matrix = matrix;
            
             % Make a Stack Plot of the spectra:

                    [datasets(d).dataStruct.plotInds,datasets(d).dataStruct.plotIndsCat] = calc_stackPlotInds({datasets(d).dataStruct.smoothedData.data},plotRes,maxTime);
                    
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
                         legend 'off'
                         
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
                                    num2str(plotRegion(1)),'-',num2str(plotRegion(2)),...
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
    




