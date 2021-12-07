%% ---- Select Data ---- %%
% This function is to select the data we are going to use in the
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
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%   2. infoListPath:The path to the data info list. (string)
%   3. originalPath: The path to the original scans. (string)
%   4. bfcPath: The path to the the bias field corrected scans. (string)
%   5. normalizedPath: The path to the normalized scans. (string)
%   6. mSegmentationPath: The path to the manual segmentations. (string)
%
% Output:
%   1. data: Tables with info about the selected cases. (struct)
%

function data = selectData(masterPath,infoListPath,originalPath,bfcPath,normalizedPath,mSegmentationPath)
%% Read the list
list = readtable(infoListPath);
%% All MRGB patients with manual segmentations
data.All = list(contains(list.PartOfMRGBSubset,'True') &...
    contains(list.ManualSegmentationAvilable,'True'),6:11);
% Remove cells with "RS" or "RE" scans (recurence) or "ST" (staiging) or
% CaseXXX (This scan has no T2W sequance)
% Loop over columns
for ii = 1:5
    columnN = data.All.Properties.VariableNames{ii+1};
    % Loop over rows
    for jj = 1:size(data.All,1)
        data.All.(columnN)(contains(data.All.(columnN),{'RE','RS','ST'}))= {''};
    end
end
% Remove cases which one of their scans have no T2W sequance
noT2W = {'Case....'};
data.All(contains(data.All.PatientID,noT2W),:) = {''};
%% Cases with multiple scans
data.MultipleScans = data.All(~cellfun(@isempty,data.All.Scan2ID),:);
%% Scanning dates
data.Date = table;
% Loop over the columns
for ii = 1:5
    columnN = data.MultipleScans.Properties.VariableNames{ii+1};
    % Loop over the rows
    for jj = 1:size(data.MultipleScans,1)
        % Convert to datetime format
        if ~isempty(data.MultipleScans.(columnN){jj})
            data.Date.(columnN)(jj) = datetime(datestr(datenum(num2str(...
                data.MultipleScans.(columnN){jj}(end-5:end),'%d'),...
                'ddmmyy'),'dd/mm/yyyy'),'InputFormat','dd/MM/yyyy');
        else
            data.Date.(columnN)(jj) = NaT('Format','defaultdate');
        end
    end
end
%% Difference between dates
data.Difference = table;
% Loop over the columns
for ii = 1:4
    columnN1 = data.Date.Properties.VariableNames{ii};
    columnN2 = data.Date.Properties.VariableNames{ii+1};
    % Measure the days difference
    data.Difference.(columnN1) = datenum(data.Date.(columnN2))-...
        datenum(data.Date.(columnN1));
end
%% Cases with difference less than 91 days
[r,c] = find(table2array(data.Difference)<91);
% Exclude repeated patients, take the first pare of scans
[~,~,idx] = unique(r,'rows');
idy = 1:numel(r);
Z = accumarray(idx(:),idy(:),[],@(n){n});
rep = find(cellfun(@numel,Z)>1);
exclude = zeros(numel(rep));
for ii = 1:numel(rep)
    exclude(ii,:) = Z{rep(ii)}(2:end);
end
exclude = unique(exclude(:));
r(exclude) = [];
c(exclude) = [];
%% Selected data
data.Selected = table;
for ii = 1:numel(r)
    data.Selected.Scan1{ii} = [char(data.MultipleScans{r(ii),c(1)}) '_' char(data.MultipleScans{r(ii),c(ii)+1})];
    data.Selected.Scan2{ii} = [char(data.MultipleScans{r(ii),c(1)}) '_' char(data.MultipleScans{r(ii),c(ii)+2})];
    data.Selected.Interval{ii} = table2array(data.Difference(r(ii),c(ii)));
end
%% Copy the selected Scans
% Output paths
scan1OF = fullfile(masterPath,'Data','Original','Scan1');
scan2OF = fullfile(masterPath,'Data','Original','Scan2');
scan1BF = fullfile(masterPath,'Data','N4BFC','Scan1');
scan2BF = fullfile(masterPath,'Data','N4BFC','Scan2');
scan1NF = fullfile(masterPath,'Data','Normalized','Scan1');
scan2NF = fullfile(masterPath,'Data','Normalized','Scan2');
scan1MSF = fullfile(masterPath,'Data','Segmentation','Manual','Original','Scan1');
scan2MSF = fullfile(masterPath,'Data','Segmentation','Manual','Original','Scan2');
% If the output files not exist make them
if ~exist(scan1OF,'dir')
    mkdir(scan1OF)
end
if ~exist(scan2OF,'dir')
    mkdir(scan2OF)
end
if ~exist(scan1BF,'dir')
    mkdir(scan1BF)
end
if ~exist(scan2BF,'dir')
    mkdir(scan2BF)
end
if ~exist(scan1NF,'dir')
    mkdir(scan1NF)
end
if ~exist(scan2NF,'dir')
    mkdir(scan2NF)
end
if ~exist(scan1MSF,'dir')
    mkdir(scan1MSF)
end
if ~exist(scan2MSF,'dir')
    mkdir(scan2MSF)
end
% Loop over the selected scans
for ii = 1:size(data.Selected,1)
    % Copy the Original scans
    copyfile(fullfile(originalPath,[data.Selected.Scan1{ii} '.mhd']),scan1OF)
    copyfile(fullfile(originalPath,[data.Selected.Scan1{ii} '.raw']),scan1OF)
    copyfile(fullfile(originalPath,[data.Selected.Scan2{ii} '.mhd']),scan2OF)
    copyfile(fullfile(originalPath,[data.Selected.Scan2{ii} '.raw']),scan2OF)
    % Copy the N4BFC scans
    copyfile(fullfile(bfcPath,[data.Selected.Scan1{ii} '_n4bfc.mhd']),scan1BF)
    copyfile(fullfile(bfcPath,[data.Selected.Scan1{ii} '_n4bfc.raw']),scan1BF)
    copyfile(fullfile(bfcPath,[data.Selected.Scan2{ii} '_n4bfc.mhd']),scan2BF)
    copyfile(fullfile(bfcPath,[data.Selected.Scan2{ii} '_n4bfc.raw']),scan2BF)
    % Copy the normalized scans
    copyfile(fullfile(normalizedPath,[data.Selected.Scan1{ii} '_Normalized.mhd']),scan1NF)
    copyfile(fullfile(normalizedPath,[data.Selected.Scan1{ii} '_Normalized.raw']),scan1NF)
    copyfile(fullfile(normalizedPath,[data.Selected.Scan2{ii} '_Normalized.mhd']),scan2NF)
    copyfile(fullfile(normalizedPath,[data.Selected.Scan2{ii} '_Normalized.raw']),scan2NF)
    % Copy the manual segmentations
    copyfile(fullfile(mSegmentationPath,data.Selected.Scan1{ii}(1:8),...
        [data.Selected.Scan1{ii} '_Segmentation.mhd']),scan1MSF)
    copyfile(fullfile(mSegmentationPath,data.Selected.Scan1{ii}(1:8),...
        [data.Selected.Scan1{ii} '_Segmentation.raw']),scan1MSF)
    copyfile(fullfile(mSegmentationPath,data.Selected.Scan2{ii}(1:8),...
        [data.Selected.Scan2{ii} '_Segmentation.mhd']),scan2MSF)
    copyfile(fullfile(mSegmentationPath,data.Selected.Scan2{ii}(1:8),...
        [data.Selected.Scan2{ii} '_Segmentation.raw']),scan2MSF)
end

%% Training CNN data
data.TrainCNN = data.All;
for ii = 1:size(data.Selected,1)
    tID = data.Selected.Scan1{ii}(1:8);
    data.TrainCNN(contains(data.All.PatientID,tID),:) = {''};
end
data.TrainCNN = data.TrainCNN(~cellfun(@isempty,data.TrainCNN.PatientID),:);
% Take the first scan to represent each patient with ine scan only
for ii = 1:size(data.TrainCNN,1)
    data.TrainCNNselected{ii,:} = [char(data.TrainCNN.PatientID(ii)) '_' char(data.TrainCNN.Scan1ID(ii))];
end
%% Copy the CNN training data scans
% Output paths
scanOF = fullfile(masterPath,'Data','Train','Original');
scanBF = fullfile(masterPath,'Data','Train','N4BFC');
scanNF = fullfile(masterPath,'Data','Train','Normalized');
scanMSF = fullfile(masterPath,'Data','Train','Segmentation','Manual','Original');
% If the output files not exist make them
if ~exist(scanOF,'dir')
    mkdir(scanOF)
end
if ~exist(scanBF,'dir')
    mkdir(scanBF)
end
if ~exist(scanNF,'dir')
    mkdir(scanNF)
end
if ~exist(scanMSF,'dir')
    mkdir(scanMSF)
end
% Loop over the selected scans
for ii = 1:size(data.TrainCNNselected,1)
    % Copy the Original scans
    copyfile(fullfile(originalPath,[data.TrainCNNselected{ii} '.mhd']),scanOF)
    copyfile(fullfile(originalPath,[data.TrainCNNselected{ii} '.raw']),scanOF)
    % Copy the N4BFC scans
    copyfile(fullfile(bfcPath,[data.TrainCNNselected{ii} '_n4bfc.mhd']),scanBF)
    copyfile(fullfile(bfcPath,[data.TrainCNNselected{ii} '_n4bfc.raw']),scanBF)
    % Copy the normalized scans
    copyfile(fullfile(normalizedPath,[data.TrainCNNselected{ii} '_Normalized.mhd']),scanNF)
    copyfile(fullfile(normalizedPath,[data.TrainCNNselected{ii} '_Normalized.raw']),scanNF)
    % Copy the manual segmentations
    copyfile(fullfile(mSegmentationPath,data.TrainCNNselected{ii}(1:8),...
        [data.TrainCNNselected{ii} '_Segmentation.mhd']),scanMSF)
    copyfile(fullfile(mSegmentationPath,data.TrainCNNselected{ii}(1:8),...
        [data.TrainCNNselected{ii} '_Segmentation.raw']),scanMSF)
end
save('data.mat','data')
end