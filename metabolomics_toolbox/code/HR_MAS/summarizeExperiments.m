function [exps,expSumm] = summarizeExperiments(exps)

%     for i = 1:length(exps)
%         fprintf([num2str(i),' - ',[exps{i}.plotTitle{:}],'\n'])
%         fprintf(['\t\t',num2str(exps{i}.density),' mg/rotor\n'])    
%         fprintf(['\t\t',num2str(length(exps{i}.smoothedData)),' samples\n'])  
%         fprintf(['\t\tRun ID for sample ',num2str(length(exps{i}.smoothedData)),' - ',exps{i}.studyInfo.sample(exps{i}.traceMat).paramFiles.name,'\n'])
%         fprintf(['\t\t',num2str(    size(exps{i}.smoothedData(exps{i}.traceMat).data,   1)  ),' timepoints in sample ',num2str(length(exps{i}.smoothedData)),'\n'])  
%     end

    for i=1:length(exps)
        expSumm(i).description = [exps{i}.plotTitle{:}];
        expSumm(i).density = exps{i}.density;
        expSumm(i).traceMat = exps{i}.traceMat;
        [~,expSumm(i).name] = fileparts(exps{i}.studyInfo.sampleDir);
        clear dts dte expnames
        for j = 1:length(exps{i}.studyInfo.sample)
            dts(j) = datetime(exps{i}.studyInfo.sample(j).paramFiles.fileData(1).date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
            dte(j) = datetime(exps{i}.studyInfo.sample(j).paramFiles.fileData(end).date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
            expnames{j} = regexprep(exps{i}.studyInfo.sample(j).paramFiles.name,[expSumm(i).name,'_'],'');
        end
        [dts,inds] = sort(dts);
        expSumm(i).startDate = string(min(dts));
        expSumm(i).endDate = string(max(dte(inds)));
        expSumm(i).experiments = strjoin(expnames(inds),', ');
        expSumm(i).smoothedSpectraPerExperiment = num2str(cellfun(@(x) size(x,1),{exps{i}.smoothedData.data}));
        exps{i}.summary = expSumm(i);
    end
    
    
end