function [sorter_settings] = get_sortersettings(filenames)

for c = 1:length(filenames);
    filereads{c,1} = fileread(filenames{c,1});
end

% using general expressions, it finds the settings at the bottom of the file
% if another setting is important to extract, a regular expression can be
% written for it and added to the list

for i = 1:length(filereads);
    pattern_flowcell = '(FlowCell)\:.(\d\d\d)';
    flowcell(i,1) = regexp(filereads{i,1},pattern_flowcell,'tokens');
    flowcell_data{i,1} = [flowcell{i,1}{1,1}];
    flowcell_data{i,2} = [flowcell{i,1}{1,2}];
    
    pattern_scanrate = '(Scan Rate)\:.(\d+)';
    scanrate(i,1) = regexp(filereads{i,1},pattern_scanrate,'tokens');
    scanrate_data{i,1} = [scanrate{i,1}{1,1}];
    scanrate_data{i,2} = [scanrate{i,1}{1,2}];
    
    pattern_triggersource = '(Trigger source)\:.(\w+)';
    triggersource(i,1) = regexp(filereads{i,1},pattern_triggersource,'tokens');
    triggersource_data{i,1} = [triggersource{i,1}{1,1}];
    triggersource_data{i,2} = [triggersource{i,1}{1,2}];
    
    pattern_signalthreshold = '(Signal threshold)\:.(\d+)';
    signalthreshold(i,1) = regexp(filereads{i,1},pattern_signalthreshold,'tokens');
    signalthreshold_data{i,1} = [signalthreshold{i,1}{1,1}];
    signalthreshold_data{i,2} = [signalthreshold{i,1}{1,2}];
    
    pattern_mintimeofflight = '(Minimum time of flight)\:.(\d+)';
    mintimeofflight(i,1) = regexp(filereads{i,1},pattern_mintimeofflight,'tokens');
    mintimeofflight_data{i,1} = [mintimeofflight{i,1}{1,1}];
    mintimeofflight_data{i,2} = [mintimeofflight{i,1}{1,2}];

    pattern_signalgains = '(Signal gains)\:\n.(Extinction)\:\s(\d.+)\n.(Green)\:\s(\d.+)\n.(Yellow)\:\s(\d.+)\n.(Red)\:\s(\d.+)\nE';
    signalgains(i,1) = regexp(filereads{i,1},pattern_signalgains,'tokens');
    signalgains_data{i,1} = [signalgains{i,1}{1,1}];
    signalgains_data{i,2} = [signalgains{i,1}{1,2}];
    signalgains_data{i,3} = [signalgains{i,1}{1,3}];
    signalgains_data{i,4} = [signalgains{i,1}{1,4}];
    signalgains_data{i,5} = [signalgains{i,1}{1,5}];
    signalgains_data{i,6} = [signalgains{i,1}{1,6}];
    signalgains_data{i,7} = [signalgains{i,1}{1,7}];
    signalgains_data{i,8} = [signalgains{i,1}{1,8}];
    signalgains_data{i,9} = [signalgains{i,1}{1,9}];
    
    pattern_ext_detector_power = '(Extinction detector power)\:.\n(\d.+\n)PMT';
    ext_detector_power(i,1) = regexp(filereads{i,1},pattern_ext_detector_power,'tokens');
    ext_detector_power_data{i,1} = [ext_detector_power{i,1}{1,1}];
    ext_detector_power_data{i,2} = [ext_detector_power{i,1}{1,2}];
    
    pattern_pmt_voltage = '(PMT voltage)\:\n.(Green)\:\s(\d.+)\n.(Yellow)\:\s(\d.+)\n.(Red)\:\s(\d\d\d)';
    pmt_voltage(i,1) = regexp(filereads{i,1},pattern_pmt_voltage,'tokens');
    pmt_voltage_data{i,1} = [pmt_voltage{i,1}{1,1}];
    pmt_voltage_data{i,2} = [pmt_voltage{i,1}{1,2}];
    pmt_voltage_data{i,3} = [pmt_voltage{i,1}{1,3}];
    pmt_voltage_data{i,4} = [pmt_voltage{i,1}{1,4}];
    pmt_voltage_data{i,5} = [pmt_voltage{i,1}{1,5}];
    pmt_voltage_data{i,6} = [pmt_voltage{i,1}{1,6}];
    pmt_voltage_data{i,7} = [pmt_voltage{i,1}{1,7}];
    
    pattern_pressure = '(Pressure)\:\n(\d....)\;(\d....)\;(\d....)';
    pressure(i,1) = regexp(filereads{i,1},pattern_pressure,'tokens');
    pressure_data{i,1} = [pressure{i,1}{1,1}];
    pressure_data{i,2} = [pressure{i,1}{1,2}];
    pressure_data{i,3} = [pressure{i,1}{1,3}];
    pressure_data{i,4} = [pressure{i,1}{1,4}];
   
    pattern_dropwidth = '(Drop width)\:.(\d\.\d)';
    dropwidth(i,1) = regexp(filereads{i,1},pattern_dropwidth,'tokens');
    dropwidth_data{i,1} = [dropwidth{i,1}{1,1}];
    dropwidth_data{i,2} = [dropwidth{i,1}{1,2}];
    
    pattern_sortdelay = '(Sort delay)\:.(\d.+)\nL';
    sortdelay(i,1) = regexp(filereads{i,1},pattern_sortdelay,'tokens');
    sortdelay_data{i,1} = [sortdelay{i,1}{1,1}];
    sortdelay_data{i,2} = [sortdelay{i,1}{1,2}];
  
end

sorter_settings = [flowcell_data scanrate_data triggersource_data signalthreshold_data mintimeofflight_data signalgains_data ...
    ext_detector_power_data pmt_voltage_data pressure_data dropwidth_data sortdelay_data];



%     pattern_gates = 'R.+\(.+\;';
%     gates(i,1) = regexp(filereads{i,1},pattern_gates,'tokens');
%     gates_data{i,1} = [gates{i,1}{1,1}];
