%% ---- Post Processing ---- %%
% This function is to post procdss the automatically segmented cases to use in the
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
% 17.Nov.2020
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%   2. quality: The quality and the calculated performance scores of the segmentations. (structure)
%
function  postProcessing(masterPath,quality)
esPath = fullfile(masterPath,'Data','Segmentation','Automated');
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
            % Prepare output path
            outP = fullfile(masterPath,'Data','Segmentation','AutomatedPP',...
                aD(ii).name,nD(jj).name,rD(kk).name);
            if ~exist(outP,'dir')
                mkdir(outP)
            end
            % Loop over cases
            sD = dir(fullfile(rD(kk).folder,rD(kk).name,'*.mhd'));
            for ll = 1:numel(sD)
                % Path in
                Fpath = fullfile(sD(ll).folder,sD(ll).name);
                % Path out
                Opath = fullfile(outP,sD(ll).name);
                % Get the largest connected component
                getLargestComp(Fpath,Opath);
            end
        end
    end
end
end

%% getLargestComp
% This function is to get the largest connected component in 3D volume and
% save it as .mhd
%
% Input:
%   1. Fpath: The full path of the segmentation mask. (string)
%   2. Opath: The full path of the output mask. (string)
%
function getLargestComp(Fpath,Opath)
% Read the scan in matlab
[StrDatax, ~, ~] = elxMetaIOFileToStrDatax(Fpath, 0);
% Find connected components in binary image
conncomp = bwconncomp(StrDatax.Data, 26);
% Identify the largest component using cellfun
[~, maxcell] = max(cellfun(@numel, conncomp.PixelIdxList));
% Zero the image and assign to it the largest component
StrDatax.Data = zeros(size(StrDatax.Data));
StrDatax.Data(conncomp.PixelIdxList{1, maxcell}) = 1;
% Save the mask as .mhd
[~, ~] = elxStrDataxToMetaIOFile(StrDatax, Opath, 0);
end