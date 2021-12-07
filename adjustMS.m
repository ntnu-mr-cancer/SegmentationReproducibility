%% ---- Adjust Manual Segmentation ---- %%
% This function is to move adjusted manual segmentations for the
% study "The Reproducibility of deep learning-based segmentation of the
% prostate on T2-weighted MR images".
%
% Sunoqrot, M.R.S.; Selnæs, K.M.; Sandsmark, E.; Langørgen, S.; Bertilsson,
% H.; Bathen, T.F.; Elschot, M. The Reproducibility of Deep Learning-Based
% Segmentation of the Prostate Gland and Zones on T2-Weighted MR Images.
% Diagnostics 2021, 11, 1690. https://doi.org/10.3390/diagnostics11091690
% https://www.mdpi.com/2075-4418/11/9/1690
%
% By Mohammed R. S. Sunoqrot, MR cancer group, NTNU, Trondheim, Norway
% 05.Oct.2020
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%
function adjustMS(masterPath)
%% The adjusted manual segmentations path
amSegmentationPath = fullfile(masterPath,'Adjusted');
%% Paths
trainP = fullfile(masterPath,'Data','Train','Segmentation','Manual','Adjusted');
scan1P =  fullfile(masterPath,'Data','SegmentationManual','Adjusted','Scan1');
scan2P =  fullfile(masterPath,'Data','Segmentation','Manual','Adjusted','Scan1');
%% Move files
movefile(fullfile(amSegmentationPath,'Train'),trainP)
movefile(fullfile(amSegmentationPath,'Scan1'),scan1P)
movefile(fullfile(amSegmentationPath,'Scan2'),scan2P)
end