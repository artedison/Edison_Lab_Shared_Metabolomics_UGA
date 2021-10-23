function [shortSpec] = reformatOldShorts(shortSpec,density,traceMat,plotRes)

                shortSpec.colorsSmoothed = shortSpec.colorsSmoothed;

                shortSpec.studyInfo = shortSpec.studyInfo;

                s = strsplit(shortSpec.studyInfo.sample(1).paths.sample,'/');
                shortSpec.studyInfo.sampleDir = s{end-1}; clear s;
                
            
                shortSpec.density = density;
                shortSpec.traceMat = traceMat;
                
                if isfield(shortSpec,'Xsmoothed')
                    if isfield(shortSpec.Xsmoothed,'data')
                        for s = 1:length(shortSpec.Xsmoothed)
                            shortSpec.smoothedData(s).data = shortSpec.Xsmoothed(s).data;
                        end
                    end
                    if isfield(shortSpec.Xsmoothed,'timepoints')
                        for s = 1:length(shortSpec.Xsmoothed)
                            shortSpec.smoothedData(s).timepoints = shortSpec.Xsmoothed(s).timepoints;
                        end
                    end
                else
                    
                    if isfield(shortSpec,'smoothedData')
                        if isfield(shortSpec.smoothedData,'data')
                            [shortSpec.trackingInds,...
                                shortSpec.trackingIndsCat,...
                                shortSpec.TrackingMatrix] = calc_stackPlotInds(shortSpec.smoothedData(shortSpec.traceMat).data,...
                                                                                plotRes,'smooth');
                        end
                        if isfield(shortSpec.smoothedData,'timepoints')
                            shortSpec.smoothedTimes(shortSpec.traceMat).timepoints = shortSpec.smoothedData(shortSpec.traceMat).timepoints;
                        end
                    end
                    
                end
                
                if ~isfield(shortSpec,'ppm')
                    shortSpec.ppm = shortSpec.ppmcat;
                end
                
            if isfield(shortSpec,'plotTitle')
                if iscell(shortSpec.plotTitle)
                    shortSpec.plotTitle = {[shortSpec.plotTitle{:}]};
                else
                    shortSpec.plotTitle = {[shortSpec.plotTitle]};
                end
            else
                shortSpec.plotTitle = {shortSpec.studyInfo.sample(1).paramFiles.name};
            end
      
                    
end