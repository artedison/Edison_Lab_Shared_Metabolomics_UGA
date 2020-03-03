%% Before you begin, make certain that the Edison Lab MATLAB metabolomics toolbox and its subfolders is added to your path. 
    % You can obtain the current version of this for free from https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
    % We recommend that you put the toolbox in some location on your local
    % computer that github syncs the toolbox too. Add all the folders and
    % subfolders to your path from that github directory.
    
%% Get the data:
%     You will need to get the processed NMR data from this location:
%     'Dropbox (Edison_Lab@UGA)/Resources/code_review/1d_nmr/MATLAB_analysis'

%% Define path
    % We recommend keeping your workflow (this file) in the directory that you
    % want to work from in a MATLAB session. It is convenient to also put your
    % workspace.mat file in the same director. 
    % The simplest way to start a session is to just double click the
    % workflow.m file, which will launch MATLAB in the right place.
    % For this tutorial, you can also read in the accompanying file "workspace_worms.mat", 
    % which has already gone through many of the steps below.
%% Generate GroupID from the NMR sample_info sheet
    % The Y vector is the way that MATLAB knows which sample is which in a study.
    % 
    % We recommend storing the sample Y vector data and other naming or metadata in a
    % folder called "NMR_info" that is one directory up from the workflow.m file. 
    % For a MS study, this folder would be called "MS_info". 
    %
    % This standardizes the directory structure from study to study, and from the
    % sample information directory, we then put the raw data in another folder on the
    % same level so that it can also be found using another file.
    %
    % If you follow these guidelines, you do not need any modification to the
    % "readtable" command below. If you chose a different option, you will need
    % to modify the filename to include the correct path.
    %
    % Read in an excel table containing information about the samples.

    % The table must have 5 columns: 
    %   Col1: Run_ID - The run ID of the sample. This should be sorted
    %                   numerically in ascending order.
    %   Col2: Sample_ID - A unique identifier for each sample.
    %   Col3: Sample_desc - This is for user's purposes. Add any description
    %                       that might be relevant to the sample for future reference.
    %   Col4: Sample_grp - The group a sample belongs to (For ex: treatment1)
    %   Col5: GroupID - Each Sample_grp gets a corresponding GroupID number. So, all
    %               samples from the same group must have the same GroupID number.
    
    SampleInfo = readtable('../NMR_info/NMR_sample_info.xlsx'); 
    GroupID = SampleInfo.GroupID;
    
%% %%%%%%%%%%%%%%% NMR DATA POST PROCESSING %%%%%%%%%%%%%%%%%%%

%% Load ft files (This command loads the processed data(.ft files from NMRPipe) into MATLAB)
    % This path also needs to be changed if the NMR_info is not one directory
    % up from the workflow.m file.

    loadStudyFTdata('../NMR_info/1D_noesypr1d_path.csv');
    
%% Setup spectra (This creates vectors that MATLAB can use to process the spectra. X is the name for the starting data vector.
    
    [X,ppm,Xtitle]=Setup1D(spectra);

%% Look at spectra to make sure they look okay. 
    % If there are big outliers that represent technical problems, they should
    % be removed here to avoid causing problems with other steps below.

    displaypeak1D(X,ppm,0,GroupID);

%% With this command you can select and label a specific spectrum
    % This helps to figure out which one has problems. Clicking a spectrum
    % displays its row number. 

    whichLine()
    
%% Referencing to DSS (you have two options doing it manually or automaticly instead of doing them one by one)
    % DSS (4,4-dimethyl-4-silapentane-1-sulfonic acid) is added to each NMR
    % sample as a reference. Here it is being used as a chemical shift
    % reference and is set to 0.00 ppm. As long as DSS is the same across all
    % samples and that it does not bind to compounds in the sample, it can also
    % be used as a quantification reference. It is important to know that the
    % signal at 0.00 ppm is comprised of 9 equivalent protons. Thus, if you add
    % 1/9 mM DSS to your sample, the integrated value of the peak at 0.00 ppm
    % will be 1 mM.
    % The second command below that automatically references DSS works well if
    % the DSS is similar in every spectrum.
    %spectra = ref_spectra_manual(spectra, [-0.02, 0.01], 0)
    
    spectra = ref_spectra(spectra,-0.02);

%% IMPORTANT NOTE: After refrencing you have to repeat the setup step

    [X,ppm,Xtitle]=Setup1D(spectra);
    
%% Check spectra to make sure all are referenced correctly
% You should zoom into the DSS region to verify that it is correctly
% referenced to 0.00 ppm and is similar in intensity and has good and
% symmetric linewidths.

    displaypeak1D(X,ppm,0,GroupID);

%% Removing ends,solvent regions and any other regions you don't want 
% This step is a "blunt hammer". Essentially, any empty or misbehaving regions can be replaced by zeros. 
% This should be done cautiously, but it is important to prevent highly variable regions from dominating the analysis.
% The settings below remove ends, water and methanol regions. They need to
% be adjusted for each sample, so never trust that these settings apply to
% generic studies.
% For each step, notice that we are updating the data vector X to new
% vectors that have had parts removed. The final vector after this step is
% XR3. 
% IMPORTANT: If you decide to skip this step, you will need to rename
% your X vector to XR3 (XR3=X) before the Alignment step below.
% Alternatively, you can change the XR3 in the "align1D" command to X.

    [XR1,ppmR] = remove_ends(X,ppm,-0.12,10.0);
    XR2=remove_region(XR1,ppmR,4.71,4.8);
    XR3=remove_region(XR2,ppmR,3.37,3.3);

%% This will show the resulting spectra with the areas defined above removed.

    displaypeak1D(XR3,ppmR,0,GroupID);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALIGNMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alignment is a "messy" step, is never perfect, and sometimes hard to even
% make it good. The reason NMR data need to be aligned is that some of the
% peak positions (chemical shift values) are quite sensitive to pH, salts,
% and other external factors. It is best to start with buffered samples,
% but that often is not enough.
%
% We offer several different alignment options through this toolbox: 
% Wong J.W.H., Durante C., Cartwright H.M. Application of fast Fourier
% transform cross-correlationfor the alignment of large chromatographic 
% and spectral datasets. Analytical Chemistry 2005;77:5655-5661.
%
% Nielsen N.P.V., Carstensen J.M., Smedsgaard J. (1998) Aligning of single 
% and multiple wavelength chromatographic profiles for chemometric data 
% analysis using correlation optimised warping. J. Chromatogr. A. 805, 17-35.
% 
%% The alignment (PAFFT) using 'correlation' seems to perform better than the other one listed below
% However, this one is not perfect. You should always explore options and
% then visualize results of different algorithms.
% 
% Notice that file after alignment is 'XAL' (X aligned). The other options
% below use XALn (n=1,9). If you want to use one of these others in the
% workflow below this step, it will need to be renamed XAL. For example, if
% 'CCOW' with 'spearman' works the best, you will need to run the command
% 'XAL=XAL6' to get the right inputfile for normalization below. 'star'
% alignment lines spectra up with a reference (often the max or mean
% spectrum). Guide align leverages HCA of the spectra, for 
% datasets where presence/absence differences are expected between groups.
%
    XAL=guide_align1D(XR3,ppmR,'spearman','PAFFT');
%
%% These other options can be tested by uncommenting the lines that start with 'XALn=...'
%
% %% Alignment (PAFFT)
% XAL1=star_align1D(XR3,ppmR,'mean','PAFFT');
% 
% %% Alignment (RAFFT)
% XAL2=star_align1D(XR3,ppmR,'mean','RAFFT');
% 
% %% Alignment (ccow)
% XAL3=star_align1D(XR3,ppmR,'mean','CCOW');
% 
% %% Alignment (ccow)
% XAL4=star_align1D(XR3,ppmR,'median','CCOW');
% 
% %% Alignment (ccow)
% XAL5=guide_align1D(XR3,ppmR,'correlation','CCOW');
% 
% %% Alignment (ccow)
% XAL6=guide_align1D(XR3,ppmR,'spearman','CCOW');
% 
% %% Alignment (RAFFT)
% XAL7=guide_align1D(XR3,ppmR,'correlation','RAFFT');
% 
% %% Alignment (RAFFT)
% XAL8=guide_align1D(XR3,ppmR,'spearman','RAFFT');
%
% %% Alignment (RAFFT)
%XAL9=guide_align1D(XR3,ppmR,'spearman','PAFFT');

%% Next you can visualize these to see which is best. 
%If you save figures (name.fig), you can compare them as a group as shown
%below with "MultiplePlot"
displaypeak1D(XAL,ppmR,0,GroupID);
%
%%
%displaypeak1D(XAL2,ppmR,0,GroupID);
%%
%displaypeak1D(XAL4,ppmR,0,GroupID);
%%
%displaypeak1D(XAL6,ppmR,0,GroupID);
%%
%displaypeak1D(XAL7,ppmR,0,GroupID);
%%
%displaypeak1D(XAL8,ppmR,0,GroupID);
%%
%displaypeak1D(XAL9,ppmR,0,GroupID);
%% Plot them together to compare and select the best alignment 
% If you have several figures (e.g. from multiple different alignments) to compare in MATLAB, you can use the following command.
% MultiplePlot is very versatile and will expand all of the figues the
% same, allowing for easy comparison. The first arguement (2 in this case)
% specifies the number of columns. 
%
% We made a folder called 'figures' that has all of the different
% alignments already made. The first 'cd' will change directory to those
% figures and the second changes back to the current directory.
%
cd 'figures/'
MultiplePlot(2,'XAL.fig','XAL1.fig','XAL2.fig','XAL3.fig','XAL4.fig','XAL5.fig','XAL6.fig','XAL7.fig','XAL8.fig','XAL9.fig');
cd '..'
%
%% Normalization
% Samples are not always the same quantity. For example, some worm samples
% could have fewer or smaller worms than others. Or in a study of human
% urine, people who drink a lot of water will have more dilute urine than
% those who don't. 
%
% Normalization attempts to make samples across a study comparable to each
% other. It is a "row operation" in that it applies the same correction to
% every point in a spectrum. It does not change the relative intensities of
% peaks within a sample but does between samples.
%
% Normcheck will show how normal the data are. A well-normalize set of
% samples will have median values that are the same.
%
%% This checks the normalization. 
%
normcheck(XAL)
%
%% Normalization by PQN 
% PQN:'The approach of probabilistic quotient normalization assumes
% that changes in the concentrations of single analytes only influence 
% parts of the spectra, whereas changes of the overall concentration of a 
% sample influence the complete spectrum. Although there are alternatives
% (e.g. see below), PQN tends to work best for most NMR metabolomics data
% sets.
%
% Dieterle, F.; Ross, A.; Schlotterbeck, G.; Senn, H. (2006) Probabilistic
% Quotient Normalization as Robust Method to Account for Dilution of
% Complex Biological Mixtures. Application to 1H NMR Metabolomics.
% Analytical Chemistry 78, 4281?4290.
%
% Note that the normalized file is called 'XALN'
%
XALN=normalize(XAL,ppmR,'PQN');
%
%% Check normalization again.
% All the data should have mean values near zero and appear uniform.
% Compare this with the previous one using 'XAL'
%
normcheck(XALN)
%
%%
whichLine()
%% Look at all the aligned and normalized data
%
displaypeak1D(XALN,ppmR,0,GroupID);
%
%% Other method for normalizations in the toolbox are 'total' for Total Area, 
%'intensity' is for normalization to single feature,'integral'
% is for normalization to sum of set of features

%% Using DSS peak to normalize the data set 
%
%XALN_DSS=normalize(XAL6,ppmR,'integral',[-0.06,0.06]);
%
%%
%displaypeak1D(XALN_DSS,ppmR,0,GroupID);
%% Plots to look at the means
% Sometimes it is most useful to just examine plots of the mean spectra.
% The set of commands below does that. It also adds linewidth (here 1). A
% nice display is to show means with a heavier linewidth and individual
% replicates with narrower linewidths in the same collor. 
%
figure
hold
plotr(ppmR,mean(XALN(GroupID==8,:)),'b','linewidth',1)
%
plotr(ppmR,mean(XALN(GroupID==7,:)),'g','linewidth',1)
%
plotr(ppmR,mean(XALN(GroupID==6,:)),'r','linewidth',1)
%
plotr(ppmR,mean(XALN(GroupID==5,:)),'m','linewidth',1)
%
plotr(ppmR,(XALN(GroupID==4,:)),'k','linewidth',1)
%
%% Filter and select subset of samples 
% There is great flexibility in MATLAB and with a Y vector to select data
% to plot or analyze. 
%
%GroupID ID
% 4 = L1
% 5 = L3
% 6 = L4
% 7 = Adults 
% 8 = Adults + hatched L1 larvae
%
selection=[6 7 8];
%
%% This plots whatever subset you want. 
% For example, try switching 'sel29' for 'sel21'
%
plot_spec(XALN,ppmR,SampleInfo,selection,'GroupID','yes','show','mean','color','mean');
% 
%% %%%%%%%%%%%%Analysis%%%%%%%%%%%
%
%% Scaling the data
% Scaling helps when large peaks dominate the analysis, and small peaks get lost. 
% An analogy is analyzing a forest for significant features. The oak trees
% will dominate, but the significant features may be the mushrooms. Scaling
% helps recognize these.
%
% Whereas normalization is a row operation, scaling is a column operation.
% Thus, it changes the relative intensity of features across a specific
% sample.
%
% There are also several options for scaling, and we include 2 of the most
% common ones here. We generally use 'logoff', but 'pareto' is often good
% too. You need to experiment to see what works best.
% 
% Here is scaling using 'logoff'. Note that the name of the matrix becomes
% XALNS, which indicates that the matrix has been aligned, normalized and
% scaled.
%
XALNS=scale(XALN,'logoff');
%
%% Scale using 'pareto'
%XALNS_P=scale(XALNFinal,'pareto');
%% Check scaling for logoff
% We generally focus on the final plot showing the distribution of features
% across points of the spectra. This should be as flat as possible.
%
varcheck(XALNS)
%
%% Principal Component Analysis
% PCA is best done with the aligned, normalized, and scaled data. You can
% explore to see what happens if you don't do that.
% This nipalsPCA command will get the first 5 components.
%
PCA=nipalsPCA(XALNS,5);
%
%% This command makes the scores plot. 
% Change the dimensions to see the PCA components of interest. 3 numbers
% (e.g. [1 2 3]) will make a 3D plot. 2 numbers will make an x-y plot.
%
figure
hold
VisScores(XALNS,PCA,[1 2],'Y',GroupID,'conf_ellipse',true);
%
%% Loadings plots: PC-1
% These show what data cause the separation along a particular PCA
% component. Notice that the data are aligned and normalized but not
% scaled. The scaling makes the loadings plots difficult to interpret.
%
VisLoadings1D(XALN,PCA.loadings(1,:),ppmR)
%
%% PC-2
%
VisLoadings1D(XALN,PCA.loadings(2,:),ppmR)
%
%% PLS-DA UNDER CONSTRUCTION.
% 7 is for 7 fold cross validation,10 isthe  number of permutations to test cross-validation
PLS=plsCV(GroupID,XALNS,7,'da',10)

%% 
VisScoresLabel(XALNS,PLS,[1 2],GroupID);
%%
VisLoadings1D(XALN,PLS.loadings(1,:),ppmR)
%%
VisLoadings1D(XALN,PLS.loadings(2,:),ppmR)
%% STOCSY
% Statistical Total Correlation Spectroscopy (STOCSY) is a powerful way to
% explore NMR data to find correlations between peaks. The strongest
% correlations (>0.9) tend to be peaks in the same molecule. One strategy
% is to take peaks from loadings in PCA or PLS-DA and STOCSY these to see
% what they are correlated with. These can then be searched in a database.
STOCSY(7.099,XALN,ppmR);
%% Test pushes