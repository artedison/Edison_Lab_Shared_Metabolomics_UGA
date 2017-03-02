%% SECTION TITLE
% DESCRIPTIVE TEXT
function [spectra, offsetppm] = ref_spectra2D(spectra,type,varargin)
% ************
% Reference spectra.
%% Required arguments: This argument is required for hte function to work
%     spectra - Structure array containing fields: real, ppm1, ppm2, Titles
%
%       example: spectra_referenced=ref_spectra2D(spectra);
%
%% Optional arguments: 
%     type - 'auto' (default) other option is 'manual'. Manual will let you
%     manually pick the referencing peaks for all spectra as opposed to
%     auto
%
%     example: spectra_referenced=ref_spectra2D(spectra,'manual'); 
%           or spectra_referenced=ref_spectra2D(spectra,'auto');
%
%% Optional Name/Value Pairs:
%     ppm= should be a vector containing ppm1 and ppm2 respectively - [0.0, 0.0] (default for both dimensions) This is the ppm that you want the shift you
%     choose to go to ok!
%
%     example: spectra_referenced=ref_spectra2D(spectra, 'ppm', [3.5 69.0]); 
%
%     ppmrange - [-1, 1; -5, 5] (default) This is the range that you will be
%     able to see your peak for peak picking
%
%     example: spectra_referenced=ref_spectra2D(spectra, 'ppmrange', [3.0 4.0; 65.0 72.0]); 
%
%% Example for a 1H - 13C 2D experiment
% [spectra_ref,offsetppm]=ref_spectra2D(spectra,'manual','ppm',[3.5 69.0],'ppmrange', [3.0 4.0; 65.0 72.0]);
%  
% Written by Chaevien Clendinen
% University of Florida, February 2015
% *******************

%% Parse Arguments
p = inputParser;
addRequired(p,'spectra',@(x)validateattributes(x,{'struct'}, {'nonempty'}));

addOptional(p, 'type', 'auto', @(x) any(validatestring(x,{'auto','manual'})));
% if exist('type','var')~=1
%     type='auto';
% end

addParamValue(p, 'ppm', [0.0, 0.0], @(x)validateattributes(x,{'numeric'}, {'vector'}));
addParamValue(p, 'ppmrange', [-1, 1; -5, 5], @(x)validateattributes(x,{'numeric'}, {'size',[2,2]}));

parse(p,spectra,type,varargin{:})

cellfun(@(f) evalin('caller',[f ' = p.Results.' f ';']), fieldnames(p.Results))

data=squeeze(spectra(1).real);
dim=dimension(data);
if dim ~= 2
    errordlg('ERROR: Script is only for 2D data! For 1D data, use ref_spectra');
    return
end

%% Manual referencing for all spectra
switch type
    case 'manual'
        disp('Manual Referencing')
        for i=1:length(spectra);
            %imagescr(spectra(i).ppm1,spectra(i).ppm2,spectra(i).real);
            vis2D(spectra(i),'f')
            xlabel('F2 (ppm1)', 'fontsize',14);
            ylabel('F1 (ppm2)', 'fontsize',14);
            title(spectra(i).Title,'fontsize',14);
            set(gca,'xlim',ppmrange(1,:),'ylim',ppmrange(2,:))
            [x,y]=ginput(1);
            
            close
            
            [~,idx(1)]=min(abs(x-spectra(i).ppm1));
            [~,idx(2)]=min(abs(y-spectra(i).ppm2));
            
            ref_ppm(i,1)=spectra(i).ppm1(idx(1));
            ref_ppm(i,2)=spectra(i).ppm2(idx(2));
            
            offsetppm(i,1)=ppm(1)-ref_ppm(i,1);
            offsetppm(i,2)=ppm(2)-ref_ppm(i,2);
            
            spectra(i).ppm1=spectra(i).ppm1+offsetppm(i,1);
            spectra(i).ppm2=spectra(i).ppm2+offsetppm(i,2);
            
        end
    case 'auto'
        %% Automated referencing for all spectra
        %Manually pick reference peak
        disp('Automated Referencing')
        i=1; %choose first spectrum
%         imagescr(spectra(i).ppm1,spectra(i).ppm2,spectra(i).real); %plot that spectrum
        vis2D(spectra(i),'f')
        title(spectra(i).Title,'fontsize',14);
        xlabel('F2 (ppm1)', 'fontsize',14);
        ylabel('F1 (ppm2)', 'fontsize',14);
        set(gca,'xlim',ppmrange(1,:),'ylim',ppmrange(2,:))
        [x,y]=ginput(1);
        close
        [~,idx(1)]=min(abs(x-spectra(i).ppm1));
        [~,idx(2)]=min(abs(y-spectra(i).ppm2));
        
        %% Automated referencing for the remainder N spectra
        if ppm(1)==0 && ppm(2)==0
            
            for i=1:length(spectra)
                
                offsetppm(i,1)=spectra(i).ppm1(idx(1));
                offsetppm(i,2)=spectra(i).ppm2(idx(2));
                
                spectra(i).ppm1=spectra(i).ppm1-offsetppm(i,1);
                spectra(i).ppm2=spectra(i).ppm2-offsetppm(i,2);
                
            end
            
        else
            
            for i=1:length(spectra)
                
                
                ref_ppm(i,1)=spectra(i).ppm1((idx(1)));
                ref_ppm(i,2)=spectra(i).ppm2((idx(2)));
                
                offsetppm(i,1)=ppm(1)-ref_ppm(i,1);
                offsetppm(i,2)=ppm(2)-ref_ppm(i,2);
                
                spectra(i).ppm1=spectra(i).ppm1-offsetppm(i,1);
                spectra(i).ppm2=spectra(i).ppm2-offsetppm(i,2);
                
            end
            
        end
        %% Plot 1D representation of spectra
        
        %For the first dimension
        subplot(1,2,1)
        xlabel('F2 (ppm1)', 'fontsize',14);
        hold on
        cmap=colormap(lines(length(spectra)));
        
        for i=1:length(spectra)
            oneDrep(i,:)=squeeze(sum(abs(spectra(i).real),1));
            plotr(spectra(i).ppm1,oneDrep(i,:),'Color',cmap(i,:))
        end
        set(gca,'xlim',[min(spectra(1).ppm1) max(spectra(1).ppm1)])
        hold off
        
        %For the second dimension
        subplot(1,2,2)
        xlabel('F1 (ppm2)', 'fontsize',14);
        hold on
        cmap=colormap(lines(length(spectra)));
        
        for i=1:length(spectra)
            oneDrep2(i,:)=squeeze(sum(abs(spectra(i).real),2));
            plotr(spectra(i).ppm2,oneDrep2(i,:),'Color',cmap(i,:))
        end
        set(gca,'xlim',[min(spectra(1).ppm2) max(spectra(1).ppm2)])
        warndlg('Displaying 1D representation of 2D referenced plot')
end
end
