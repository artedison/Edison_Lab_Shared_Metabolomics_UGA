function [fig,ax] = plotSpectra(X,ppm,varargin)

    % MTJ 2020
    titleStr = '';
    xlab = 'Chemical Shift δ (ppm)';
    ylab = 'Spectral Intensity';
    
    if ~isempty(varargin)
        
        [~,ind] = ismember('title',varargin);
        if any(ind)
            titleStr = varargin{find(ind,1)+1};
        end
        
        [~,ind] = ismember('xlabel',varargin);
        if any(ind)
            xlab = varargin{find(ind,1)+1};
        end
        
        [~,ind] = ismember('ylabel',varargin);
        if any(ind)
            ylab = varargin{find(ind,1)+1};
        end
        
    end
    
    fig = figure;
    ax = plotr(ppm,X);
    xlabel(xlab)
    ylabel(ylab)
    title(titleStr)
end