
regexp(fdata,['(?<=nmrPipe -in )\S*','\s\\'],'match')

regexp(fdata,['(?<!#[^\n])',... % not preceded by a # followed by any non-newline char (e.g. not commented out)
                '(?<=nmrPipe -in )',... % preceded by 'nmrPipe -in '
                '\S*',...   % expression is any number of consecutive non-whitespace chars
                '(?=\s\\)'],... % followed by ' \'
                'match') 

regexp(fdata,['(?<!\n#\s*)-out'],'match')
regexp(fdata,['(?<!\n#\s*)(?<=-out )','\S*','(?= -ov)'],'match')
regexp(fdata,['(?<!#[^\n])(?<=-out )','\S*','(?= -ov)'],'match') % not preceded by a # followed by any non-newline char (e.g. not commented out)

replaceIn = 'inStr';
replaceOut = 'outStr';
fdata = regexprep(fdata,['(?<!#[^\n])',...  % not preceded by a # followed by any non-newline char (e.g. not commented out)
                          '(?<=nmrPipe -in )',... % preceded by 'nmrPipe -in '
                            '\S*',...       % expr is any number of consecutive non-whitespace chars
                          '(?=\s\\)'],...   % followed by ' \'
                         replaceIn);        % replacement for expr
            
fdata = regexprep(fdata,['(?<!#[^\n])',... % not preceded by a # followed by any non-newline char (e.g. not commented out)
                          '(?<=-out )',... % preceded by '-out '
                            '\S*',...      % expr is any number of consecutive non-whitespace chars
                          '(?= -ov)'],...  % followed by ' -ov'
                         replaceOut);      % replacement for expr
