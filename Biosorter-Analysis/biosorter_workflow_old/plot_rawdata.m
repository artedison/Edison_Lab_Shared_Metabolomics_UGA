%% selects the folder the .txt files are saved in
%# build a list of file names with absolute path
fPath = uigetdir('.', 'Select directory');
if fPath==0, error('no folder selected'), 
end
fNames = dir( fullfile(fPath,'*.txt') );
fNames = strcat(fPath, filesep, {fNames.name});
%% imports and erases any cells that are empty

% preallocating array space
TOF = cell(1, length(fNames));
Ext = cell(1, length(fNames));
Green = cell(1, length(fNames));
Red = cell(1, length(fNames));
Yellow = cell(1, length(fNames));

% imports arrays
for i=1:length(fNames)
   TOF{i} = import_tof(fNames{i});
   Ext{i} = import_ext(fNames{i});
   Green{i} = import_green(fNames{i});
   Red{i} = import_red(fNames{i});
   Yellow{i} = import_yellow(fNames{i});
end

% erases empty cells
fh = @(x) all(isnan(x(:)));

TOF(cellfun(fh, TOF)) = [];
Ext(cellfun(fh, Ext)) = [];
Green(cellfun(fh, Green)) = [];
Red(cellfun(fh, Red)) = [];
Yellow(cellfun(fh, Yellow)) = [];

% erases empty cells

TOF1 = cell(1, length(TOF));
for k = 1:length(TOF);
TOF1{1,k} = num2cell(TOF{1,k});
TOF1{1,k}(cellfun(fh, TOF1{1,k})) = [];
end

Ext1 = cell(1, length(Ext));
for k = 1:length(Ext);
Ext1{1,k} = num2cell(Ext{1,k});
Ext1{1,k}(cellfun(fh, Ext1{1,k})) = [];
end

Green1 = cell(1, length(Green));
for k = 1:length(Green);
Green1{1,k} = num2cell(Green{1,k});
Green1{1,k}(cellfun(fh, Green1{1,k})) = [];
end

Red1 = cell(1, length(Red));
for k = 1:length(Red);
Red1{1,k} = num2cell(Red{1,k});
Red1{1,k}(cellfun(fh, Red1{1,k})) = [];
end
 
Yellow1 = cell(1, length(Yellow));
for k = 1:length(Yellow);
Yellow1{1,k} = num2cell(Yellow{1,k});
Yellow1{1,k}(cellfun(fh, Yellow1{1,k})) = [];
end
%%
TOF2 = cell(1, length(TOF1));
for k = 1:length(TOF1);
TOF2{1,k} = cell2mat(TOF1{1,k});
%[TOF2{1,k}] = deleteoutliers(TOF2{1,k});
end

Ext2 = cell(1, length(Ext1));
for k = 1:length(Ext1);
Ext2{1,k} = cell2mat(Ext1{1,k});
%[Ext2{1,k}] = deleteoutliers(Ext2{1,k});
end
%%
for hh = 1:59;
    Ext2_CS1{1,hh} = Ext2{1,hh};
    TOF2_CS1{1,hh} = TOF2{1,hh};
    subplot (10,6,hh);
    plot(TOF2_CS1{1,hh},Ext2_CS1{1,hh},'bo','Color','blue')
    xlim([0 1600]);
    ylim([0 2000]);
end
