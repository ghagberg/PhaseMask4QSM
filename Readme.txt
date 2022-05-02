        %============== BMMR University Tuebingen ===============%
        %========= High Field MRI Max-Planck-Institute ==========%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  Supplementary DATA                    %  
        %   how to generate masks including brain surface voxels %
        %                  described in                          %
        %'Phase-based masking for QSM of the human brain at 9.4T'%
        %Hagberg GE et al, MagnResonMed 2022                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The provided material can be used in different ways

1. Code: QSM2016challenge_MM_PB_OrigMasks.m
This code can be run in a Linux environment provided that you have installed FSL and set your paths to the following toolboxes (please change the relevant lines in the .m code for this):
addpath(genpath('PathToSPM12GoesHere'));
addpath(genpath('PathToSTI_SuiteGoesHere'));
addpath(genpath('PathToMEDI_toolboxGoesHere'));

It first takes the magnitude image (./FromChallenge/magn_raw.nii)  and generates a magnitude-based mask (MM) with bet
then the wrapped phase ('./FromChallenge/phs_wrap.nii') from the QSM2016challenge is taken to identify PB 
These two approaches, proposed in the paper, can be compared with the provided mask ('./FromChallenge/msk.nii')

after laplacian unwrapping and V-SHARP background removal, qsm images are generated using different dipole algorithms

N.B substitute the magn_raw and phs_wrap with your own (correctly scaled) images and try out the code on your own GRE data.


2. if you do not want to run the code you can just browse the  Results folders.

ResultsLapl 
ResultsRomeo 

They contain:

a. unwrapped phase images after background removal: 
uwRSH*.nii using Laplacian and RESHARP Tk-12
roRSH*.nii using ROMEO and RESHARP Tk-12
uwVSH*.nii using Laplacian and V-SHARP, SMV20mm
roVSH*.nii using ROMEO and V-SHARP, SMV20mm

N.B. ROMEO uses a mask as input, so it was run three times, to compare performance of ROMEO unwrapping using differ3n masks.

b. qsm images using different dipole inversion algorithms:
qsmL2*  using the L2-algorithm provided in the challenge
qsmLSQR* from STI-Suite
qsmMEDI* from the MEDI toolbox

3. Performance metrics for each method have been extracted from a mask valid for all techniques (as reported in Supplementary Table S2) and are available in TableSupplDat.xls
The code FinalTests4SuppDat.m describes how this was done

4. Code: ExtractROIdat_kFe.m
This code used to determine the iron-dependent QSM-contrast. For this purpose, the QSM2016 challenge data has been segmented into CSF, grey, and white matter. 
There is also a subdivision into cortical and subcortical areas according to the Desikan-Killiany (Harvard-Oxford atlas in FSL) atlas.
The tissue and ROI segmentated images are available in the folder: kFE
Based on the reported age of the subject, region-specific estimates of non-haeme iron concentration in microgram Iron per gram wet weight tissue is provided in the code.
Linear regression, and the statistics for the regression curve is obtained The QSM-contrast, kFe, is the slope of the curve 
