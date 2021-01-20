function [response,doExit,doCancel] = CIVM_processingMenu()

    % Present the user with a menu of options

        response = menu('Pick one of the following options:"',... % 
                    'Show data',...                               % 1
                    'Phase Spectra',...                           % 2
                    '',...                          % 3
                    '',...                      % 4
                    'Apply Processing to Dataset',...             % 5
                    'Reset phasing (in development)',...                           % 6
                    'Restart from ft.com template upload',...     % 7
                    'SAVE AND EXIT',...                           % 8
                    '',...
                    'Cancel (do not save)');                      % 10
                
    % Inform on Save and exit options
    
        doExit = false;
        doCancel = false; 
        
        if response == 8
            doExit = true;
        elseif response == 10
            doCancel = true;
        end
end