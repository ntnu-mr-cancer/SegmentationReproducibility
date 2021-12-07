%% ---- Prepare automated segmentation masks ---- %%
% This function is to prepare the automated segmentations masks
% that we are going to use in the
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
% 23.Nov.2020
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%
function prepAS(masterPath)
%% Paths
% The main output folder
asPath = fullfile(masterPath,'Data','Segmentation','Automated');
if ~exist(asPath,'dir')
    mkdir(asPath)
end
% The automated segmentations path
aSegmentationPath = fullfile(masterPath,'Automated');
%% Types
type1 = {'V-Net'};
type2 = {'nnU-Net'};
%% Type 1
%-- Copy and correct WP and PZ masks
% Loop over the segmenation networks
aD = dir(aSegmentationPath);
aD = aD(~ismember({aD.name},{'.','..'}));
aD = aD(~contains({aD.name},type2));
for ii = 1:numel(aD)
    % Network output folder
    oN = fullfile(asPath,aD(ii).name);
    if ~exist(oN,'dir')
        mkdir(oN)
    end
    % Loop over the region
    nD = dir(fullfile(aSegmentationPath,aD(ii).name));
    nD = nD(~ismember({nD.name},{'.','..'}));
    for jj = 1:numel(nD)
        % Region output folder
        oR = fullfile(oN,nD(jj).name);
        if ~exist(oR,'dir')
            mkdir(oR)
        end
        % Loop over scans
        rD = dir(fullfile(aSegmentationPath,aD(ii).name,nD(jj).name));
        rD = rD(~ismember({rD.name},{'.','..'}));
        for kk = 1:numel(rD)
            % Scan output folder
            oS = fullfile(oN,nD(jj).name,rD(kk).name);
            if ~exist(oS,'dir')
                mkdir(oS)
            end
            % Loop over cases
            sD = dir(fullfile(aSegmentationPath,aD(ii).name,nD(jj).name,rD(kk).name));
            sD = sD(~ismember({sD.name},{'.','..'}));
            sD = sD(~contains({sD.name},'tested'));
            for ll = 1:numel(sD)
                source = fullfile(aSegmentationPath,aD(ii).name,...
                    nD(jj).name,rD(kk).name,sD(ll).name);
                target = fullfile(oS,sD(ll).name);
                if ~isfile(target)
                    % to make it readable with matlab later and save in
                    % target
                    [StrData, ~, ~] = elxMetaIOFileToStrDatax(source, 0);
                    [~, ~] = elxStrDataxToMetaIOFile(StrData, target, 0);
                end
            end
        end
    end
end
%-- Create TZ
% Loop over the segmenation networks
for ii = 1:numel(aD)
    % Loop over the scans
    scan = {'Scan1','Scan2'};
    for jj = 1:numel(scan)
        % Loop over cases
        wpD = dir(fullfile(asPath,aD(ii).name,'WP',scan{jj},'*.mhd'));
        pzD = dir(fullfile(asPath,aD(ii).name,'PZ',scan{jj},'*.mhd'));
        tzD = fullfile(asPath,aD(ii).name,'TZ',scan{jj});
        if ~exist(tzD,'dir')
            mkdir(tzD)
        end
        for kk = 1:numel(wpD)
            wpF = fullfile(wpD(kk).folder,wpD(kk).name);
            pzF = fullfile(pzD(kk).folder,pzD(kk).name);
            tzF = fullfile(tzD,wpD(kk).name);
            if ~isfile(tzF)
                % Read WP and PZ
                [wp, ~, ~] = elxMetaIOFileToStrDatax(wpF, 0);
                [pz, ~, ~] = elxMetaIOFileToStrDatax(pzF, 0);
                tz = wp;
                tz.Data(pz.Data==1) = 0;
                [~, ~] = elxStrDataxToMetaIOFile(tz, tzF, 0);
            end
        end
    end
end

%% Type 2
% Loop over the segmenation networks
aD = dir(aSegmentationPath);
aD = aD(~ismember({aD.name},{'.','..'}));
aD = aD(~contains({aD.name},type1));
for ii = 1:numel(aD)
    % Network output folder
    oN = fullfile(asPath,aD(ii).name);
    if ~exist(oN,'dir')
        mkdir(oN)
    end
    % Loop over scans
    rD = dir(fullfile(aSegmentationPath,aD(ii).name));
    rD = rD(~ismember({rD.name},{'.','..'}));
    for kk = 1:numel(rD)
        % Scan output folder
        oS = fullfile(masterPath,'Temp',aD(ii).name,rD(kk).name);
        mkdir(oS)
        % Target output paths
        oSwp = fullfile(oN,'WP',rD(kk).name);
        oSpz = fullfile(oN,'PZ',rD(kk).name);
        oStz = fullfile(oN,'TZ',rD(kk).name);
        if ~exist(oSwp,'dir')
            mkdir(oSwp)
        end
        if ~exist(oSpz,'dir')
            mkdir(oSpz)
        end
        if ~exist(oStz,'dir')
            mkdir(oStz)
        end
        % Activate environment to use python
        if isunix
            setenv('PATH','.../Anaconda/anaconda3/envs/qcs/bin');
        end
        % Wrtie paths to txt file to read in python
        fileID = fopen(fullfile(masterPath,'paths.txt'),'w');
        fprintf(fileID,'%s\n',...
            fullfile(aSegmentationPath,aD(ii).name,rD(kk).name));
        fprintf(fileID,'%s\n',oS);
        fclose(fileID);
        % Use SimpleITK from python to conver nii.gz to mhd and
        % copy file
        [~,~] = system("python niigztomhdnnUNet.py");
        % reset default environment
        if isunix
            setenv('PATH','/usr/bin');
            rehash toolboxcache;
            savepath
        end
        delete paths.txt
        % Correct masks header, Delete some tags to enable ElastixFromMatlab toolbox
        % Loop over the cases
        dirlist = dir(fullfile(oS,'*.mhd'));
        A = cell(1,1);
        for mm = 1:numel(dirlist)
            filename = fullfile(dirlist(mm).folder,dirlist(mm).name);
            % Read mhd header into cell A
            fid = fopen(filename,'r');
            i = 1;
            tline = fgetl(fid);
            A{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fid);
                A{i} = tline;
            end
            fclose(fid);
            % Change cell A and save in B
            B = {A{1:10},A{58:60}}; % new cell array with the wanted lines
            % Write cell B into mhd header
            fid = fopen(filename, 'w');
            for i = 1:numel(B)
                fprintf(fid,'%s\n', B{i});
            end
            fclose(fid);
            % Get WP, PZ and TZ masks and save them to there target
            [StrDatax, ~, ~] = elxMetaIOFileToStrDatax(filename, 0);
            WP = StrDatax;
            WP.Data = zeros(size(WP.Data));
            PZ = StrDatax;
            PZ.Data = zeros(size(PZ.Data));
            TZ = StrDatax;
            TZ.Data = zeros(size(TZ.Data));
            % Get masks
            WP.Data(StrDatax.Data>0) = 1;
            PZ.Data(StrDatax.Data==1) = 1;
            TZ.Data(StrDatax.Data==2) = 1;
            % Save
            [~, ~] = elxStrDataxToMetaIOFile(WP, fullfile(oSwp,dirlist(mm).name), 0);
            [~, ~] = elxStrDataxToMetaIOFile(PZ, fullfile(oSpz,dirlist(mm).name), 0);
            [~, ~] = elxStrDataxToMetaIOFile(TZ, fullfile(oStz,dirlist(mm).name), 0);
        end
    end
end
% Remove Temp folder
rmdir(fullfile(masterPath,'Temp'),'s')
end