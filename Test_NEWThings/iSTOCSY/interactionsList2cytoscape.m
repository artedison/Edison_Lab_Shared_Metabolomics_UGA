function [] = interactionsList2cytoscape(interactions,net_outFile)

%% Reformat interactions for Cytoscape (binary interactions)
% Takes in: 
% interactions: list of interactions
% exclude: list of ppms to exclude
% net_outFile: name of output file containing binary interactions

[rows,~] = size(interactions);
% for i = 1:row, row = take the the data from the row without zeros or NAN
  for i = 1:rows;
      [~,~,row] = find(interactions(i,:));
    for j = 2:length(row), 
        %print [row(1),rowlength (j);
        dlmwrite(net_outFile, [interactions(i,1),interactions(i,j)],'-append');
    end
  end

interactions = csvread(net_outFile);

% Remove unwanted ppm values from the interactions
%{
for i = 1:length(exclude)
    interactions = interactions(~or((interactions(:,1)==9.168),(interactions(:,2)==9.168)),:);
    interactions = interactions(~or((interactions(:,1)==9.4383),(interactions(:,2)==9.4383)),:);
end
%}
csvwrite(net_outFile,interactions);
%Read in that which is not noise (9.168 and 9.4383 in this case) before moving on (use Cytoscape, circle
%view, get list of 'noise' (one node with a big circle around it). Just
%want a vector of numbers.

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

