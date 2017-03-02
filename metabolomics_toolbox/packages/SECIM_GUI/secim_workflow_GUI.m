%% loading a template of parameter table
table_of_parameters=readtable('chosen_workflow_params.csv');
% setting up our parameters table
parameters=parameter(table_of_parameters);
%% loading all ft files based on their path
loadStudyFTdata('path_for_ft_files.csv');
%% finding region for referencing GUI
reference_spectra_GUI(parameters, spectra); 
%% referencing spectra
spectra = ref_spectra_manual(spectra,[parameters.T.ref_thresh_1,parameters.T.ref_thresh_2]);
%% 
[X,ppm,Xtitles] = Setup1D(spectra);
%% region removal GUI
remove_region_GUI(parameters, X, ppm); 
%% removing chosen regions
X_waterRemoved = remove_region(X,ppm,parameters.T.water_region_1,parameters.T.water_region_2);
[X_waterAndEndsRemoved, ppmEndsRemoved]=remove_ends(X_waterRemoved,ppm,parameters.T.ends_1,parameters.T.ends_2);
X_RegionsRemoved=remove_region(X_waterAndEndsRemoved,ppmEndsRemoved,parameters.T.region_removed_1,parameters.T.region_removed_2); 
dlmwrite('referenced_baselineCorrected_regionRemoved.txt',X_RegionsRemoved)
dlmwrite('ppmOriginal.txt',ppm);
dlmwrite('ppmRemoved.txt',ppmEndsRemoved);
X=X_RegionsRemoved;
ppm=ppmEndsRemoved;
%% baseline correction GUI
baseline_GUI(parameters,ppm,X); 
%% correcting baseline with chosen parameters
X_baseline_corrected = CorrectBl(X, str2num(parameters.T.bl_noise), str2num(parameters.T.bl_fitting));
dlmwrite('referenced_baselineCorrected.txt',X_baseline_corrected )
X=X_baseline_corrected;
%% Alignment GUI
alignment_GUI(parameters,ppm,X);
%% Aligning with chosen parameters
if isequal(parameters.T.align_function,'Star')
    X_aligned=star_align1D(X,ppm,parameters.T.align_represent,parameters.T.align_method);
else
    X_aligned=guide_align1D(X,ppm,parameters.T.align_represent,parameters.T.align_method);
end
dlmwrite('aligned.txt', X_aligned);
X=X_aligned;
%% Normalization GUI
Normalization_GUI(parameters,ppm,X);
%% Normalizing with the chosen parameters
if isequal(parameters.T.normalization_feature,'NA')
    X_normalized=normalize(X,ppm,parameters.T.normalization_method);
else
    X_normalized=normalize(X,ppm,parameters.T.normalization_method,parameters.T.normalization_feature);
end
dlmwrite('Normalized.txt', X_normalized);
X=X_normalized;
%% Scaling GUI
scale_GUI(parameters,X);
%% Scaling with the chosen method
X_scaled=scale(X,parameters.T.scaling_method);
dlmwrite('Normalized.txt', X_scaled);
X=X_scaled;
%% writing chosen parameters into a csv file
writetable(parameters.T,'~/matlab-test/NMR/chosen_workflow_params.csv');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%        EXTRA TESTING FOLLOWS       %%%%%%%%%%%%%%%%%%%
