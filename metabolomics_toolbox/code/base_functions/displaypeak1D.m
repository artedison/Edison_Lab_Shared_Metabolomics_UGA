function displaypeak1D(X,ppm,shift,Y)

    % Author: Edison Lab
    % Ver 0.3
    % Tested on Matlab Version R2017b
    % Date: 4FEB2019
    %
    %
    % Description:
    %   displaypeak plots a set of chemical shift values as circles on a set of
    %   1D spectra, where the colors of the spectra correspond to GroupID.
    %   Use this script to assign correlated peaks.
    %
    % Input:
    %       X: stack of 1D spectra
    %       ppm: chemical shift vector
    %       shift: one or more chemical shifts, takes scalar or vector 
    %       Y: Optional- if response vector Y is known, use it to color the scores plot
    %
    % Output:
    %       A plot of the 1D spectra.
    %
    % Log:
    %   Edited by : Rahil Taujale
    %   Date      : 2MAR2017
    %   Ver       : 0.2
    %           Changed colormap to generate it using distinguishable_colors
    %           function.
    %   Ver       : 0.3
    %           Fixed bug (Get correct vector size for colors based on Y)
    %           MJ 24JAN2018 map values in Yvec to avoid indexing problems with
    %           Y(k) = 0.
    % Example run:
    %

    % Check whether group ID is given
        if exist('Y')==0
            Y=ones(size(X,1),1);
        end
    
    % Build colormap
        [map,~,ind] = unique(Y); % MJ 24JAN2018 map values in Yvec to avoid indexing problems with Y(k) = 0.
        cmap=distinguishable_colors(length(map));
        
    % Generate the figure
        figure, hold on;        
        for k=1:size(X,1)
            plot(ppm,X(k,:),'Color',cmap(ind(k),:))
        end
        set(gca,'XDir','reverse')
        for k=1:length(shift)
            [a,b(k)]=min(abs(ppm-shift(k)));
        end
        scatter(shift,max(X(:,b)));
    
end