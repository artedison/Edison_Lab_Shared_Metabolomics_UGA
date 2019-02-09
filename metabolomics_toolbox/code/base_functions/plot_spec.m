function plot_spec(X,ppm,T,sel,varargin)

% Author: Rahil Taujale
% Version: 0.5
% Date: 16 NOV 2018
%
% Description:
%       plot_spec plots a set of chemical shift values on a set of 1D
%       spectra. It has added functionalities for display selection,
%       coloring and displaying labels.
%
% Input: 
%       X: stack of 1D spectra (n x m matrix)
%       ppm: chemical shift vector (1 x m vector)
%       T: Table with the sample information (n rows)
%           5 columns with following exact headers:
%                Run_ID,Sample_ID,Sample_desc,Sample_grp,Yvec
%       sel: 3 options:
%           1) A cell array with a list of Sample_grp names to be plotted
%               Eg: sel={'L1','L2'}
%           2) A numeric vector with a list of Yvec to be plotted
%               Eg: sel=[5 6 7]
%           3) 'all' If want to plot all samples
%       Optional:
%           'Yvec'  : 'Yes' if sel input above is a list of Yvec as in option
%                      2
%                   : 'No' (Default) if sel input is a cell array with
%                      Sample_grp names
%           'show'  : 'all' (Default) to show all the plots
%                   : 'sample' to show only the sample plots
%                   : 'mean' to show only the mean plots
%           'color' : 'mean' (Default) Label and color plots based on means
%                   : 'all' Label each sample with separate color
%
% Output: 
%       A plot of the 1D spectra (figure object/axes)
%
% Log:
%       Ver 0.2: Added secondary labels to label added peaks.
%       Ver 0.3: Removed secondary labels option. Use data_cursor_labels.m
%       function
%       Ver 0.4: When plotting, does not use the black color anymore.
%       Ver 0.5: Code Review annotated.

% Example run: plot_spec(XALN,ppmR,T,selected_sets,'Yvec','n','show','mean','color','mean')


%% Read in parameters for varargin
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
    
%% Set up the color groups depending on sel option 

    % Setup based on Yvec
        if strcmpi(sel,'all')
            sel=unique(T.Yvec)';
            opt.Yvec='yes';
        end
        
    % Generating N+1 distinguishable colors
        [~,N]=size(sel);
        CM = distinguishable_colors(N+1);   
        if (N>3)
            CM(4,:)=[];
        end
        
    % Initialize other vars
        fig=[];
        ax=[];
        sample_name=cell(0);
        figure;
        hold on;
    
    % Loop through all samples to group them
    
        for f=1:N
            % Group current sample/row based on the selection made above
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
            
            % Plotting 
            
                if (D>1)
                    [fig(end+1),ax(end+1)]=plotr(ppm,mean(Xsel),'color',CM(f,:),'linewidth',2);
                elseif (D==1)
                    [fig(end+1),ax(end+1)]=plotr(ppm,Xsel,'color',CM(f,:),'linewidth',2);
                end
                m(f)=length(fig);
                sample_name(end+1)=sel2(f);
                
        end
        
%% Adding labels, colors and legend

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

   