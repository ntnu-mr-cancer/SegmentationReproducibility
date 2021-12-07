%% ---- Feature Extraction ---- %%
% This function is to extract the shape features using the
% automated segmentations masks, to use in the
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
% 30.Oct.2020
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%
% Output:
%   1. features: The extracted shape features. (structure)
%
function features = featureExtraction(masterPath)
%% Go to the folders to extract features
% Loop over the segmentation netwroks
msPath = fullfile(masterPath,'Data','Segmentation','Manual','Final');
% asPath = fullfile(masterPath,'Data','Segmentation','Automated');
asPath = fullfile(masterPath,'Data','Segmentation','AutomatedPP');
aD = dir(asPath);
aD = aD(~ismember({aD.name},{'.','..'}));
aD(numel(aD)+1).name = 'Manual'; % add manual
for ii = 1:numel(aD)
    % Loop over the regions
    if strcmp(aD(ii).name,'Manual')
        nD = dir(fullfile(asPath,aD(1).name));
        nD = nD(~ismember({nD.name},{'.','..'}));
    else
        nD = dir(fullfile(asPath,aD(ii).name));
        nD = nD(~ismember({nD.name},{'.','..'}));
    end
    % Loop over the scans
    for jj = 1:numel(nD)
        % Loop over Scans
        if strcmp(aD(ii).name,'Manual')
            rD = dir(fullfile(asPath,aD(1).name,nD(jj).name));
            rD = rD(~ismember({rD.name},{'.','..'}));
        else
            rD = dir(fullfile(asPath,aD(ii).name,nD(jj).name));
            rD = rD(~ismember({rD.name},{'.','..'}));
        end
        for kk = 1:numel(rD)
            % Features path
            fPath = fullfile(masterPath,'Data','Features',...
                aD(ii).name,nD(jj).name,rD(kk).name);
            if ~exist(fPath,'dir')
                mkdir(fPath)
            end
            % Wrtie paths to txt file to read in python
            fileID = fopen(fullfile(masterPath,'paths.txt'),'w');
            % Images folder
            fprintf(fileID,'%s\n',...
                fullfile(masterPath,'Data','Normalized',rD(kk).name));
            % Masks folder
            if strcmp(aD(ii).name,'Manual')
                fprintf(fileID,'%s\n',fullfile(msPath,rD(kk).name,nD(jj).name));
            else
                fprintf(fileID,'%s\n',fullfile(asPath,aD(ii).name,nD(jj).name,rD(kk).name));
            end
            % Features output folder
            fprintf(fileID,'%s\n',fPath);
            % Region Class
            fprintf(fileID,'%s\n',nD(jj).name);
            fclose(fileID);
            %% Extract features
            % Activate 'srep' environment to use python
            if isunix
                setenv('PATH','.../Anaconda/anaconda3/envs/qcs/bin');
            end
            % Use Pyradiomics (V 3.0) package from python (3.7)
            %---Shape 3D
            % Wrtie paths to txt file to read in python
            [~,~] = system("python pyradiomicsFeaturesExtraction3D.py");
            % reset default environment
            if isunix
                setenv('PATH','/usr/bin');
                rehash toolboxcache;
                savepath
            end
            delete paths.txt
        end
    end
end
%% Organize features
% Loop over networks
fsPath = fullfile(masterPath,'Data','Features');
fD = dir(fsPath);
fD = fD(~ismember({fD.name},{'.','..'}));
for ii = 1:numel(fD)
    % Loop over the regions
    nD = dir(fullfile(fsPath,fD(ii).name));
    nD = nD(~ismember({nD.name},{'.','..'}));
    for jj = 1:numel(nD)
        % Loop over the Scans
        rD = dir(fullfile(fsPath,fD(ii).name,nD(jj).name));
        rD = rD(~ismember({rD.name},{'.','..'}));
        for kk = 1:numel(rD)
            % current directory name
            dir_name = fullfile(fsPath,fD(ii).name,nD(jj).name,rD(kk).name);
            % region classes
            region_classes = nD(jj).name;
            % feature classes
            feature_classes = {'shape'};
            % cases names
            namedir = dir(fullfile(masterPath,'Data','Original',rD(kk).name,'*.mhd'));
            name = {namedir.name}.';
            %% fill arrays
            for ff = 1:length(feature_classes)
                % create table
                table_name = sprintf('table_%s_%s',region_classes,feature_classes{ff});
                eval([table_name ' = table;']);
                for pp = 1:numel(name)
                    % get jsonname
                    json_name = sprintf('%s_%s_%s.json',name{pp, 1}(1:end-4),region_classes,feature_classes{ff});
                    % read in and convert
                    if exist(fullfile(dir_name,json_name),'file')
                        data = loadjson(fullfile(dir_name,json_name));
                        names = fieldnames(data);
                        fdata = rmfield(data,names(~contains(names,feature_classes{ff})));
                        colNames = fieldnames(fdata);
                        % preallocate cell array
                        colNamesNew = cell(numel(colNames),1);
                        % loop to get variables names array
                        for cc = 1:numel(colNames)
                            colNamesNew{cc,:} = colNames{cc}(10:end);
                        end
                        % create temp table for each case
                        tmp_table = struct2table(fdata);
                        tmp_table.Properties.RowNames = {sprintf(name{pp, 1}(1:end-4))};
                        tmp_table.Properties.VariableNames = colNamesNew;
                        % add the tables together
                        eval([table_name ' = [' table_name '; tmp_table];']);
                        clear colNames colNamesNew
                    end
                end
                % if not empty, save file
                if not(isempty(eval(table_name)))
                    % replace the Network name, replace "-" by "_"
                    netName = strrep(fD(ii).name,'-','_');
                    % save features to structure
                    features.(netName).(nD(jj).name).(rD(kk).name) = eval(table_name);
                end
            end
        end
    end
end
% Save features structure
save('features.mat','features')
end