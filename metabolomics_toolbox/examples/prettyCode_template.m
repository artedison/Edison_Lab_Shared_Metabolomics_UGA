function functionName(args)

    % Author: Rahil Taujale
    % Version: 0.4
    % Tested on Matlab Version R2017b
    % Date: 21MAR2018
    %
    % Description:
    %       plot_spec plots a set of chemical shift values on a set of 1D
    %       spectra. It has added functionalities for display selection,
    %       coloring and displaying labels.
    %
    % Input:
    %       X: stack of 1D spectra
    %       ppm: chemical shift vector
    %       T: Table with the sample information.
    %           5 columns with following exact headers:
    %                Run_ID,Sample_ID,Sample_desc,Sample_grp,Yvec
    %
    % Output:
    %       A plot of the 1D spectra.
    %
    % Log:
    %       Ver 0.4: When plotting, does not use the black color anymore.
    %
    % Example run:
    %
    
    %% This is what this section does
        opt = struct('Yvec','No','show','all','color','mean','xlabel',0);
        optNames = fieldnames(opt);
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
            error('Input Argument-Value in pairs')
        end
        for pairs = reshape(varargin,2,[])
            inpName = pairs{1};
            if any(strcmp(inpName,optNames))
                opt.(inpName) = pairs{2};
            else
                error('%s not a valid parameter',inpName)
            end
        end
    %% Do something else
        if strcmpi(sel,'all')
            sel=unique(T.Yvec)';
            opt.Yvec='yes';
        end
        [~,N]=size(sel);
        CM = distinguishable_colors(N+1);
        if (N>3)
            CM(4,:)=[];
        end
        fig=[];
        ax=[];
        sample_name=cell(0);
        figure; hold on;

        % For loop for something
            for f=1:N
                if strncmpi(opt.Yvec, 'yes',1)
                    Xsel=X(T.Run_ID(T.Yvec==sel(f)),:);
                    [D,~]=size(Xsel);
                    sample_name(end+1:end+D)=T.Sample_ID(T.Yvec==sel(f));
                    for i=1:N
                        sel2{i}=char(unique(T.Sample_grp(T.Yvec==sel(i))));
                    end
                else
                    Xsel=X(T.Run_ID(strcmp(T.Sample_grp, sel(f))),:);
                    [D,~]=size(Xsel);
                    sample_name(end+1:end+D)=T.Sample_ID(strcmp(T.Sample_grp, sel(f)));
                    sel2=sel;
                end
                [fig(end+1:end+D),ax(end+1:end+D)]=plotr(ppm,Xsel,'color',CM(f,:),'linewidth',0.5,'linestyle','--');
                if (D>1)
                    [fig(end+1),ax(end+1)]=plotr(ppm,mean(Xsel),'color',CM(f,:),'linewidth',2);
                elseif (D==1)
                    [fig(end+1),ax(end+1)]=plotr(ppm,Xsel,'color',CM(f,:),'linewidth',2);
                end
                m(f)=length(fig);
                sample_name(end+1)=sel2(f);
        end

        if strcmpi(opt.show, 'sample')
            set(ax(m),'Visible','off');
        elseif strcmpi(opt.show, 'mean')
            set(ax,'Visible','off');
            set(ax(m),'Visible','on');
            opt.label='mean';
        end

        if strcmpi(opt.color, 'all')
            CM = distinguishable_colors(length(ax)+length(m));
            for i=1:length(ax)
                j=i+length(m);
                set(ax(i),'color',CM(j,:));
            end
            for i=1:length(m)
                set(ax(m(i)),'color',CM(i,:));
            end
            legend(ax,sample_name);
        elseif strcmpi(opt.color, 'mean')
            legend(ax(m),sample_name(m));
        end
        hold off;
end

