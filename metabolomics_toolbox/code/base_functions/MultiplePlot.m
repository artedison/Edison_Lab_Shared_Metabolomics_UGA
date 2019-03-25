function MultiplePlot(varargin)

    % Author: Rahil Taujale
    % Version: 0.1
    % Date: 23AUG2016
    %
    % Description:
    %       MultiplePlot allows viewing multiple spectra plots saved as MATLAB
    %       .fig files side by side in one single plot. The axes are linked so
    %       that same zooming operation occurs in all panes for easier
    %       comparison of plots.
    %
    % Input:
    %       1st argument: A number specifying the number of columns desired in
    %                     the output.
    %       2nd-nth argument: List of .fig files to be viewed
    %
    % Output:
    %       A plot with all the fig files shown side by side.
    %
    % Log:
    %       MJ edit 9OCT2017: Line 41: add the title for each plot
    %
    % Example Run:
    %       MultiplePlot(2,'Fig1.fig','Fig2.fig');

K=varargin{1};
    N=nargin-1;
    for m=1:N
        f=m+1;
        fig(m)=hgload(varargin{f});
        ax(m)=gca;
    end
    figure;
    for i=1:N
        % create and get handle to the subplot axes
        s(i) = subplot(ceil(N/K),K,i); 
        % get handle to all the children in the figure
        copyobj(allchild(get(fig(i),'CurrentAxes')),s(i));
        allchild(fig(i));
        xlab=get(get(ax(i),'xlabel'),'string');
        ylab=get(get(ax(i),'ylabel'),'string');
        xlabel(xlab);ylabel(ylab);
        title(ax(i).Title.String)
        hold on
    end
    linkaxes(s, 'xy');
    set(s,'xdir','reverse');
end
