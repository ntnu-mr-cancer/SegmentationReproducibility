%% ---- Master ---- %%
% This script is to perform all of the required analysis for the
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
% 22.Sep.2020
%
%% Clean
close all
clear
clc
warning off
%% Add Paths
addpath(genpath('.../codes'))
addpath(genpath(pwd))
%% Master Path
masterPath = '.../projects/SegmentationReproducibility';
%% Settings
settings.DataSelection = 1;
settings.AdjustManualSegmentation = 1;
settings.PrepareNewManualmasks = 1;
settings.PrepareAutomatedSegmentations = 1;
settings.QualityControl = 1;
settings.PostProcessing = 1;
settings.FeaturesExtraction = 1;
settings.StatsticalAnalysis = 1;
settings.ReportResults = 1;
%% Data Selection
disp('Step 1: Data Selection')
if settings.DataSelection
    if ~exist('data.mat','file')
        % Paths
        % The path to the original t2w images where all of them there in
        % mhd format
        originalPath = '.../mhd';
        % The path to the images which are normalized using AutoRef method
        % https://github.com/ntnu-mr-cancer/AutoRef
        normalizedPath = '.../Normalized';
        % The path to the images those N4 bias field corected
        % Check the code in n4BiasFieldCorrection.m to apply it
        bfcPath = '.../BFCorrected';
        % The path to the folder where you have manually segmented cases
        % in mhd format
        mSegmentationPath = '.../manualSegmentation';
        % The path to an excel sheet that tells the names of scan 1 and
        % scan 2 for each of the patients
        infoListPath = '.../List.xlsx';
        % Select data
        data = selectData(masterPath,infoListPath,...
            originalPath,bfcPath,normalizedPath,mSegmentationPath);
    else
        load('data.mat')
    end
end
disp('        FINISHED')
%% Adjust Manual Segmentation
disp('Step 2: Adjust Manual Segmentation')
if settings.AdjustManualSegmentation
    adjustMS(masterPath)
end
disp('        FINISHED')
%% Prepare new manual masks
disp('Step 3: Prepare New Manual masks')
if settings.PrepareNewManualmasks
    prepNMS(masterPath);
end
disp('        FINISHED')
%% Prepare Automated Segmentations
disp('Step 4: Prepare Automated Segmentations')
if settings.PrepareAutomatedSegmentations
    prepAS(masterPath);
end
disp('        FINISHED')
%% Quality Controling
disp('Step 5: Quality Controling')
if settings.QualityControl
    load('factors_InHouse.mat')
    esPath = fullfile(masterPath,'Data','Segmentation','Automated');
    quality = qualityControl(masterPath,esPath,factors,1);
else
    load('quality.mat')
end
disp('        FINISHED')
%% Post Processing
disp('Step 6: Post Processing')
if settings.PostProcessing
    postProcessing(masterPath,quality);
end
disp('        FINISHED')
%% Features Extraction
disp('Step 7: Features Extraction')
if settings.FeaturesExtraction
    features = featureExtraction(masterPath);
else
    load('features.mat')
end
disp('        FINISHED')
%% Statstical Analysis
disp('Step 8: Statstical Analysis')
if settings.StatsticalAnalysis
    load('factors_InHouse.mat')
    esPath = fullfile(masterPath,'Data','Segmentation','AutomatedPP');
    stats = statisticalAnalysis(masterPath,esPath,features,factors);
else
    load('stats.mat')
end
disp('        FINISHED')
%% Report Results
disp('Step 9: Report Results')
if settings.ReportResults
    reportResults(masterPath,stats);
end
disp('        FINISHED')