function printCleanPDF(fig,title,varargin)
%% printCleanPDF
% I've always had issues printing things to pdf without having to crop,
% etc., and MATLAB's functions are a little unintuitive. This is a 
% collection of code from:
%   https://www.mathworks.com/matlabcentral/answers/311820-save-a-figure-as-pdf
% that will allow you to tighten up the figure and print only what you see
% in the figure window, as you see it in the figure window, as a vectorized
% pdf. 
%
% "The first two lines measure the size of your figure (in inches). 
% The next line configures the print paper size to fit the figure size. 
% The last line uses the print command and exports a vector pdf document as the output."

% fillpage = '';
% vectorized = '';
% title = num2str(now());
% 

% if ~isempty(varargin)
%     '-fillpage'
%     '-vectorized' -> '-painters'
%     'tight' -> if statement
%     'title' - get next varargin
% end

    % Decide on units
        fig = gcf;
        set(fig,'Units','inches');
    
    % Get info about how the figure appears on your screen
        screenposition = get(fig,'Position'); 
        
    % Set the export settings to exactly match what's in the figure window.
    % This removes the need to crop later. 
    
    % if 'tight'
        set(fig,...                           
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);

    % Print it out.
    
        print(fig,...
        title,...
         '-dpdf','-painters')
% 
%     % Fill page instead of cropping
%     
%         fig=gcf;
%         fig.PaperPositionMode='auto';
%         print(title,'-dpdf',fillpage)

end