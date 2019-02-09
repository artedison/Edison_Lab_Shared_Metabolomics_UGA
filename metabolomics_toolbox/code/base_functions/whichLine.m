function whichLine()
  % Author: GSS (Greg Stupp)
    % Version: 0.1
    % Tested on Matlab Version R2017b
    % Date: 4FEB2019
    %
    % Description:
    %       Activates the data cursor in the active figure window to
    %       provide the row number of the spectrum that is clicked. Depends
    %       on whichLineCallback().
    %       Note: see selectLine() to get a vector output for this
    %       function.
    %
    % Input:
    %       None required, but there should be an active figure window with data. 
    %
    % Output:
    %       None, but the row number is displayed on the figure.
    %
    % Log:
    %
    % Example run:
    %

dcm = datacursormode(gcf);
datacursormode on;
set(dcm,'UpdateFcn',@whichLineCallback);
end