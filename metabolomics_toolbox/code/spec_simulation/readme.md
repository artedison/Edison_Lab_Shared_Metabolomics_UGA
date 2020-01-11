# NMR spectra simulation based on gissmo library

## structure

lib_const.m: gissmo library download and construction

nmr_spec_simu.m: function for nmr spectra simulation

simu_script.m: the main simulation script

spec_conc_simu.m: function for simulating concentration matrix

test_script.m: test script


## data

gissmo library: Resources/gissmo_lib

->"wholegissmo_spectral.complex.more.mat": spectra

->"wholegissmo_spectral.real.list.mat": peak list

->"gissmo_tab.mat": compound information table

simulation data: Projects/Bioinformatics_modeling/spectral.related/spec_simu_gissmo (you can work in your folder rather than here)

## usage

Try to read and modify simu_script.m, especially the directory path. The main part to tune the simulation is the "USER ARGUMENT" block. If you have questions, suggestions, or bugs, email me (yue.wu@uga.edu)
