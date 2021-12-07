%% ---- Quality Controling ---- %%
% This function is to control the quality of the segmented cases to use in the
% study "The Reproducibility of deep learning-based segmentation of the
% prostate on T2-weighted MR images".
%
% Sunoqrot, M.R.S.; Selnæs, K.M.; Sandsmark, E.; Langørgen, S.; Bertilsson,
% H.; Bathen, T.F.; Elschot, M. The Reproducibility of Deep Learning-Based
% Segmentation of the Prostate Gland and Zones on T2-Weighted MR Images.
% Diagnostics 2021, 11, 1690. https://doi.org/10.3390/diagnostics11091690
% https://www.mdpi.com/2075-4418/11/9/1690
%
% This is a published work, which has a public code at
% https://github.com/ntnu-mr-cancer/SegmentationQualityControl
%
% In case of using or refering to this system, please cite it as:
% Sunoqrot, M.R.S.; Selnæs, K.M.; Sandsmark, E.; Nketiah, G.A.;
% Zavala-Romero, O.; Stoyanova, R.; Bathen, T.F.; Elschot, M.
% A Quality Control System for Automated Prostate Segmentation on
% T2-Weighted MRI. Diagnostics 2020,10, 714.
% https://doi.org/10.3390/diagnostics10090714
%
% https://www.mdpi.com/2075-4418/10/9/714
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%   2. esPath: The path to the automatically segmented cases. (string)
%   3. factors: The factors we need to transfer the mertics to percentage. (structure)
%   4. saveFlag: A logical flag to decide if you want to save the outpu. (logical)
%
% output:
%   1. quality: The quality and the calculated performance scores of the segmentations. (structure)
%
function  quality = qualityControl(masterPath,esPath,factors,saveFlag)
%% Calculate the perfromance scores
quality.pScores = calculateScores(masterPath,esPath,factors);
%% Get the quality classes
quality.classes = getClasses(quality.pScores);
%% Save
if saveFlag
    save('quality.mat','quality')
end
end

%% calculateScores
% This function is to calculate the segmentation performance scores
% for multiple casea
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%   2. esPath: The path to the automatically segmented cases. (string)
%   3. factors: The factors we need to transfer the mertics to percentage. (structure)
%
% output:
%   1. scores: The calculated performance scores for the cases. (structure)
%
function scores = calculateScores(masterPath,esPath,factors)
% Master Segmentations path
msPath = fullfile(masterPath,'Data','Segmentation','Manual','Final');
% Normalized images path
rinPath = fullfile(masterPath,'Data','Normalized');
% Loop over networks
aD = dir(esPath);
aD = aD(~ismember({aD.name},{'.','..'}));
for ii = 1:numel(aD)
    % Get regions names
    nD = dir(fullfile(esPath,aD(ii).name));
    nD = nD(~ismember({nD.name},{'.','..'}));
    for jj = 1:numel(nD)
        % Get scans names
        rD = dir(fullfile(nD(jj).folder,nD(jj).name));
        rD = rD(~ismember({rD.name},{'.','..'}));
        % Loop over scans
        for kk = 1:numel(rD)
            % Make table to fill
            % replace the Network name, replace "-" by "_"
            netName = strrep(aD(ii).name,'-','_');
            % Loop over cases
            sD = dir(fullfile(rD(kk).folder,rD(kk).name,'*.mhd'));
            for ll = 1:numel(sD)
                rmP = fullfile(msPath,rD(kk).name,nD(jj).name,sD(ll).name);
                emP = fullfile(sD(ll).folder,sD(ll).name);
                riP = fullfile(rinPath,rD(kk).name,[sD(ll).name(1:end-17) '_Normalized.mhd']);
                % Get Quality score
                [aQualityScore,~] = SQC(riP,emP,2,85);
                % Get score
                score = getScore(factors,rmP,emP);
                % Loop over regions structures inside score struct
                regs = fieldnames(score);
                regs = regs(~contains(regs,'Score') & ~contains(regs,'middle'));
                for mm = 1:numel(regs)
                    % Fields inside structure
                    fr = fieldnames(score.(regs{mm}));
                    % Loop over them
                    for nn = 1:numel(fr)
                        % Assign values
                        scores.(netName).(nD(jj).name).(rD(kk).name).([regs{mm} 'Score'])(ll,:) =...
                            score.([regs{mm} 'Score']);
                        scores.(netName).(nD(jj).name).(rD(kk).name).totalScore(ll,:) =...
                            score.totalScore;
                        scores.(netName).(nD(jj).name).(rD(kk).name).(regs{mm}).(fr{nn})(ll,:) =...
                            score.(regs{mm}).(fr{nn});
                        scores.(netName).(nD(jj).name).(rD(kk).name).names{ll,:} = sD(ll).name(1:end-17);
                        scores.(netName).(nD(jj).name).(rD(kk).name).qualityScore(ll,:) =...
                            aQualityScore;
                    end
                end
                % clear
                clear score
            end
        end
    end
end
end

%% getScore
% This function is to calculate the segmentation performance score
% for one case
%
% Input:
%   1. factors: The factors we need to transfer the mertics to percentage. (structure)
%   2. rmP: The referance mask path. (string)
%   3. emP: The estimated mask path. (string)
%
% output:
%   1. score: The calculated performance scores for one case. (structure)
%
function score = getScore(factors,rmP,emP)
% Read the masks
% MetaIO to Structured data
[referenceMaskStrData, ~, ~] = elxMetaIOFileToStrDatax(rmP, 0);
[estimatedMaskStrData, ~, ~] = elxMetaIOFileToStrDatax(emP, 0);
% Structured data to 3D image structure
referenceMask = elxStrDataxToIm3d(referenceMaskStrData);
estimatedMask = elxStrDataxToIm3d(estimatedMaskStrData);
% Region classes
region_classes = {'WP','apex','middle','base'};
% Assign slices based on region class
% Reference mask
% Get the slices contains mask
rNrVoxels = squeeze(sum(sum(referenceMask.Data)));
rfSlices = find(rNrVoxels>0);
rSlices = rfSlices(1):rfSlices(end);
% Assign the regions' slices
rNSl = length(rSlices);
rdNSl = rNSl/3;
rfNSl = fix(rdNSl);

rNS.WP = rSlices(1:length(rSlices));
rNS.apex = rSlices(1:rfNSl);
rNS.middle = rSlices(rfNSl+1:rNSl-rfNSl);
rNS.base = rSlices(rNSl-rfNSl+1:rNSl);
% Slices to be excluded when deal with a specific region
rExclude.WP = [];
rExclude.apex = [rNS.middle,rNS.base];
rExclude.middle = [rNS.apex,rNS.base];
rExclude.base = [rNS.apex,rNS.middle];
% Estimaed mask
% Get the slices contains mask
eNrVoxels = squeeze(sum(sum(estimatedMask.Data)));
efSlices = find(eNrVoxels>0);
eSlices = efSlices(1):efSlices(end);
% Assign the regions' slices
eNSl = length(eSlices);
edNSl = eNSl/3;
efNSl = fix(edNSl);

eNS.WP = eSlices(1:length(eSlices));
eNS.apex = eSlices(1:efNSl);
eNS.middle = eSlices(efNSl+1:eNSl-efNSl);
eNS.base = eSlices(eNSl-efNSl+1:eNSl);
% Slices to be excluded when deal with a specific region
eExclude.WP = [];
eExclude.apex = [eNS.middle,eNS.base];
eExclude.middle = [eNS.apex,eNS.base];
eExclude.base = [eNS.apex,eNS.middle];
% Calculate the metrics
% Reference mask
rm = referenceMask;
rmR = rm.Data;
% Esitmated mask
em = estimatedMask;
emR = em.Data;
% Loop over the regions
for ii = 1:numel(region_classes)
    region_class = region_classes{ii};
    % Reference mask input
    rm.Data = rmR;
    rm.Data(:,:,rExclude.(region_class)) = 0;
    % Esitmated mask input
    em.Data = emR;
    em.Data(:,:,eExclude.(region_class)) = 0;
    % % Feed to getMetrics function
    score.(region_class) =...
        perfMetrics(rm,em,rmR,emR,rExclude.(region_class),eExclude.(region_class));
    clear rm.Data em.Data
end
% Calculate scores from metrics
% WP
WP.DSC = max([ones(1,1),[score.WP.DSC].']*factors.wholeprostate_DSC,0);
WP.aRVD = max([ones(1,1),abs([score.WP.aRVD].')]*factors.wholeprostate_aRVD,0);
WP.HD95 = max([ones(1,1),[score.WP.HD95].']*factors.wholeprostate_HD95,0);
WP.ASD = max([ones(1,1),[score.WP.ASD].']*factors.wholeprostate_ASD,0);

score.WPScore = mean([WP.DSC,WP.aRVD,WP.HD95,WP.ASD],2);
if score.WPScore > 100
    score.WPScore = 100;
elseif score.WPScore < 0
    score.WPScore = 0;
end
% apex
apex.DSC = max([ones(1,1),[score.apex.DSC].']*factors.apex_DSC,0);
apex.aRVD = max([ones(1,1),abs([score.apex.aRVD].')]*factors.apex_aRVD,0);
apex.HD95 = max([ones(1,1),[score.apex.HD95].']*factors.apex_HD95,0);
apex.ASD = max([ones(1,1),[score.apex.ASD].']*factors.apex_ASD,0);

score.apexScore = mean([apex.DSC,apex.aRVD,apex.HD95,apex.ASD],2);
if score.apexScore > 100
    score.apexScore = 100;
elseif score.apexScore < 0
    score.apexScore = 0;
end
% base
base.DSC = max([ones(1,1),[score.base.DSC].']*factors.base_DSC,0);
base.aRVD = max([ones(1,1),abs([score.base.aRVD].')]*factors.base_aRVD,0);
base.HD95 = max([ones(1,1),[score.base.HD95].']*factors.base_HD95,0);
base.ASD = max([ones(1,1),[score.base.ASD].']*factors.base_ASD,0);

score.baseScore = mean([base.DSC,base.aRVD,base.HD95,base.ASD],2);
if score.baseScore > 100
    score.baseScore = 100;
elseif score.baseScore < 0
    score.baseScore = 0;
end
% total score
score.totalScore = mean([WP.DSC,WP.aRVD,...
    WP.HD95,WP.ASD,...
    apex.DSC,score.apex.aRVD,apex.HD95,apex.ASD,...
    base.DSC,score.base.aRVD,base.HD95,base.ASD],2);
if score.totalScore > 100
    score.totalScore = 100;
elseif score.totalScore < 0
    score.totalScore = 0;
end
end

%% perfMetrics
% This function is to calculate the segmentation performance metrics
%
% Input:
%   1. referenceMask: The referance mask in im3d format. (structure)
%   2. estimatedMask: The estimated mask in im3d format. (structure)
%   3. rmR: The referance mask data in im3d format. (double)
%   4. emR: The estimated mask data in im3d format. (double)
%   5. rEX: The excluded slices from the referance mask. (double)
%   6. eEX: The excluded slices from the estimated mask. (double)
%
% output:
%   1. pMetrics: The calculated performance metrics. (structure)
%
function pMetrics = perfMetrics(referenceMask,estimatedMask,rmR,emR,rEX,eEX)
%% Define
% Estimated mask (Automted)voxels
ES = estimatedMask.Data(:);
% Sum of the estimated mask voxels
sES = sum(ES);
% Ground truth mask voxels
GT = referenceMask.Data(:);
% Sum of the ground truth mask voxels
sGT = sum(GT);
% Indexs with voxels
rmF = find(referenceMask.Data==1);
emF = find(estimatedMask.Data==1);
% if only one of the masks not there then metric to worst
if xor(numel(rmF)==0,numel(emF)==0)
    pMetrics.DSC = 0;
    pMetrics.aRVD = 100;
    pMetrics.HD95 = 100;
    pMetrics.ASD = 100;
    % if both of the masks there then calculate the metric
else
    % Dice similarity coefficient (DSC)
    pMetrics.DSC = dice(GT,ES);
    % Absolute Relative volume difference (RVD)
    pMetrics.aRVD = abs(((sES/sGT)-1)*100); % according to Heinmann et al
    % Surface Distance metrics
    % this section is copied and modified from https://github.com/emrekavur/CHAOS-evaluation/blob/master/Matlab
    %--Extract border voxels
    % Esitmated mask
    fES = emR & ~imerode(emR,strel('sphere',1));
    fES(:,:,eEX) = 0;
    fESIdx = find(fES==1);
    [x1,y1,z1] = ind2sub(size(fES),fESIdx);
    BorderVoxelsES = [x1,y1,z1];
    % Ground truth mask
    fGT = rmR & ~imerode(rmR,strel('sphere',1));
    fGT(:,:,rEX) = 0;
    fGTIdx = find(fGT==1);
    [x2,y2,z2] = ind2sub(size(fGT),fGTIdx);
    BorderVoxelsGT = [x2,y2,z2];
    %--Transforms index points to the real world coordinates
    % Esitmated mask
    realPointsES = zeros(size(BorderVoxelsES,1),size(BorderVoxelsES,2));
    for i = 1:size(BorderVoxelsES,1)
        P = estimatedMask.A*[BorderVoxelsES(i,1),BorderVoxelsES(i,2),BorderVoxelsES(i,3),1]';
        realPointsES(i,:) = P(1:3)';
    end
    % Ground truth mask
    realPointsGT = zeros(size(BorderVoxelsGT,1),size(BorderVoxelsGT,2));
    for i = 1:size(BorderVoxelsGT,1)
        P = referenceMask.A*[BorderVoxelsGT(i,1),BorderVoxelsGT(i,2),BorderVoxelsGT(i,3),1]';
        realPointsGT(i,:) = P(1:3)';
    end
    %--Distance between border voxels
    % Esitmated mask
    MdlKDT_ES = KDTreeSearcher(realPointsES);
    [~,distIndex1] = knnsearch(MdlKDT_ES,realPointsGT);
    distIndex1 = distIndex1';
    % Ground truth mask
    MdlKDT_GT = KDTreeSearcher(realPointsGT);
    [~,distIndex2] = knnsearch(MdlKDT_GT,realPointsES);
    distIndex2 = distIndex2';
    % 95% Maximum Symmetric Surface Distance (MSSD) / Haussdorf distance (95%)
    pMetrics.HD95 = max([prctile(distIndex1,95),prctile(distIndex2,95)]);
    % Average Symmetric Surface Distance (ASD)/average boundary distance (ABD)
    pMetrics.ASD = (sum(distIndex1)+sum(distIndex2))/(size(distIndex1,2)+size(distIndex2,2));
end
end

%% getClasses
% This function is to get the quality classrs of the segmentations
%
% Input:
%   1. scores: The calculated performance scores for the cases. (structure)
%
% output:
%   1. classes: The quality classes of the cases segmentation. (structure)
%
function classes = getClasses(scores)
% Loop over Netwroks
nets = fieldnames(scores);
for ii = 1:numel(nets)
    % Loop over regions
    regs = fieldnames(scores.(nets{ii}));
    for jj = 1:numel(regs)
        % Loop over scans
        scans = fieldnames(scores.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            % Get the classes
            classes.(nets{ii}).(regs{jj}).(scans{kk}).totalScore =...
                scores.(nets{ii}).(regs{jj}).(scans{kk}).totalScore>80;
            classes.(nets{ii}).(regs{jj}).(scans{kk}).WPScore =...
                scores.(nets{ii}).(regs{jj}).(scans{kk}).WPScore>80;
%             classes.(nets{ii}).(regs{jj}).(scans{kk}).qualityScore =...
%                 scores.(nets{ii}).(regs{jj}).(scans{kk}).qualityScore>80;
        end
    end
end
end