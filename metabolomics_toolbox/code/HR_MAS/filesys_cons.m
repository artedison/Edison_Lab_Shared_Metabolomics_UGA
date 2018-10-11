function filesys_cons(goaldir,datadir,sampledir)
  %
  % filesys_cons(goaldir,datadir,sampledir)
  % this function construct the template folder for systematically arrange NMR data
  % It will copy the raw data to the folder structure
  % it will also copy useful script to the folder
  % list of data that exist
  % goaldir: the location to create the project folder
  % datadir: the location to fetch data; the folder is expected to have multiple samples which contains multiple time points
  % sampledir: the location that the sample script is from
  %%yue wu revise 10/10/2018
  
  % the function will move new folders that does not exist there.
  cd(goaldir);
  mkdir 'sampleGroup';
  cd('./sampleGroup');
  % copyfile(strcat(sampledir,'*.m'),'./');
  % copyfile(strcat(sampledir,'*.com'),'./');
  % copyfile(strcat(sampledir,'*.txt'),'./');
  dirs=dir(datadir);
  dirsname={dirs.name};
  dirsname=dirsname(~ismember(dirsname,{'.','..'}));
  existsname=dir(strcat(goaldir,'sampleGroup'));
  existsname={existsname.name};
  existsname=existsname(~ismember(existsname,{'.','..'}));
  dirsname=setdiff(dirsname,existsname);
  for i = 1:length(dirsname)
    cd(strcat(goaldir,'sampleGroup'))
    % expname=strcat(subdatadirname,num2str(i));
    expname=char(dirsname(i));
    mkdir(expname);
    cd(expname);
    expdir=pwd();
    mkdir('data');mkdir('results');mkdir('scripts');
    cd('./data');
    mkdir('nmrpipe');mkdir('raw');
    cd('./nmrpipe');
    mkdir('fid');mkdir('ft');
    cd('../raw');
    copyfile(cell2mat(strcat(datadir,dirsname(i),'/*')),'./');
    cd(expdir);
    cd('./scripts');
    copyfile(strcat(sampledir,'*.m'),'./');
    copyfile(strcat(sampledir,'*.com'),'./');
    copyfile(strcat(sampledir,'*.txt'),'./');
  end
end
