# peak grouping based on ridge tracing

## Method

H-K sign is used to find ridge point which is complemented by a selected and small amount of maximum point. These point are interesected with local maximum to screen possible ridge points. Power two orthogonal polynomial are used in the computation process.

These points are connected into local segments of close points, which are connected into ridges. The angle theta between segment vector are controled below threhold thedtheta. Repeats for same time point are eliminated by selecting the point with the largest intensity.

Interactive ridge picking are used to finally screen results.

## Usage
Pleas follow the script test.new.rid.trac.multispec.workscript.m

Even though there is six tunning parameters, including  thredseg, maxaddon, lengthrid, threhold_K, threhold_H, and windsize, only thredseg and maxaddon are recommended to be modified. The idea is fixed value for other parameters can work most of the time. thredseg is the maxmum distance to connect points into a segment. The peak that move in ppm can be distant in ppm distance, and increase this value when the peak move a lot. maxaddon is for the purpose of supplement peaks when the peak is not classified as ridge. It will add local maximum from high to low in intensity. For a give window, the high peak can be forced to be ridge. The procedure of the program can be adjusted by peakflag, remflag, and totalautoflag.

Before running any workflow, the Metabolic toolbox (WEBLINK) need to be downloaded and added to MATLAB path.

1. Run the test run to check reasonable value of thredseg, and apply that value to all samples
2. In each tracing, the user will be guided through a process:
    - Delete ridges: delete ridges if the user feel that ridges is both not good traced and will affect tracing of other peaks, especially for the refine step(next step). Ridges can be deleted by clicking "Delete Ridge(s)" in the dialogs and clicking the ridges to delete. To finish the process, the user need to type enter key.
    - Refine highly overlapped region: select the region to refine ridges. The refinement process will be assume straight for each ridges to solve mistake in ridge tracing for overlapped region. The region to refine needs to be selected by a box. If there is no need for refine =, then just click the figure. Check the result to see whether the peak is as expected. smalwid_thredseg can be increased with caution to help the ridge on tip. This step can be done multiple times.
    - Pick final ridges: clicking "Pick Final Clusters" in the dialogs and clicking the ridges (remember the order) to select. To finish the process, the user need to type enter key. In the input dialog, input name for ridges(peaks) in order you click and separated by ';'.
    - End refinement: if you click yes on the dialog, then you can select lines(end by return) and draw a box (end by click) to remove those end region that are not good traced. This process can be ended by space after a whole process (line selection and box drawing). This can be used to ensure that remaining region is clean and good traced

## Attention
1. the format of result is intended to be similar to the previous version, though there is differences becuse of different methods.
2. Overlap region is still hard to deal with. It is recommended that throw away results in clearly overlap region. However, the decision is still based on the user.
3. Choose a smaller region might help ridge tracing
4. To make sure that ridges/peaks match for different sample, space row need to be added manually.
5. The corresponding workflow can be found in the folder ridgetracing_paper_workflow/
