%% ---- Statstical Analysis ---- %%
% This function is to calculate the performance and evaluation metrics
% those extacted from Scan 1 and Scan 2 for the
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
% 02.Nov.2020
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%   2. esPath: The path to the automatically segmented cases post processing. (string)
%   3. features: The extracted shape features. (structure)
%   4. factors: The factors we need to transfer the mertics to percentage. (structure)
%
% output:
%   1. stats: The statistical analysis results. (structure)
%
function stats = statisticalAnalysis(masterPath,esPath,features,factors)
%% Statistical test type
test = {'signrank','ranksum'};
%---- Main investigation
%% Get Scanning protocol parameters for investigation set
stats.protocol = getProtocol(masterPath);
%% Get Scanning protocol parameters for training set
stats.protocolTr = getProtocolTr(masterPath);
%% Calculate ICC
stats.ICC = calculateICC(features);
%% Calculate the significancy for ICC
stats.sICC = calculateSicc(stats.ICC);
%% Calculate CV
stats.CV = calculateCV(features);
%% Calculate the significancy for CV
stats.sCV = calculateScv(stats.CV,test{1});
%% Calculate the perfromance scores
stats.pScores = qualityControl(masterPath,esPath,factors,0);
%% Calculate the significancy for performance metrics
stats.sScores = calculateSscores(stats.pScores.pScores,test{1});
%% Calculate number of included slices
stats.nSlices = calculateNslices(masterPath,esPath,[]);
%% Calculate the significancy for number of included slices
stats.sSlices = calculateSslices(stats.nSlices,test{1});
%% Calculate the volume
stats.Volume = calculateVolume(features);
%% Calculate the significancy for volume
stats.sVolume = calculateSvolume(stats.Volume,test{1});
%% Calculate the change in volume
stats.cVolume = calculateCvolume(features);
%% Calculate significancy forthe change in volume
stats.sCvolume = calculateSCvolume(stats.cVolume,test{1});
%% Calculate Median +/- SD for Scores
stats.MSscores = clalcMSscores(stats.pScores.pScores);
%% Calculate Median +/- SD for number of slices
stats.MSslices = clalcMSslices(stats.nSlices);
%% Calculate Median +/- SD for Change in Volume
stats.MScvolume = clalcMScvolume(stats.cVolume);
%% Calculate DSC Difference Median +/- SD
stats.DMSscores = clalcDMSscores(stats.pScores.pScores);

%---- After Implementing QC system
%% Get Idx of the excluded cases based on the quality control score
stats.ExcludedIdx = excludedIdx(stats.pScores.pScores);
%% Get features structure after excluding cases
stats.FeaturesAE = featuresAE(features,stats.ExcludedIdx);
%% Calculate ICC
stats.ICCAE = calculateICC(stats.FeaturesAE);
%% Generate ICC after 1000 random exclusion
stats.ICCR = iccRN(features,stats.ExcludedIdx);
%% Calculate the significancy for ICC
stats.sICCAE = calculateSicc(stats.ICCAE);
%% Calculate CV
stats.CVAE = calculateCV(stats.FeaturesAE);
%% Calculate the significancy for CV
stats.sCVAE = calculateScv(stats.CVAE,test{2});
%% Calculate the perfromance scores
stats.pScoresAE = getScoresAE(stats.pScores.pScores,stats.ExcludedIdx);
%% Calculate the significancy for performance metrics
stats.sScoresAE = calculateSscores(stats.pScoresAE,test{2});
%% Calculate number of included slices
stats.nSlicesAE = calculateNslices(masterPath,esPath,stats.ExcludedIdx);
%% Calculate the significancy for number of included slices
stats.sSlicesAE = calculateSslices(stats.nSlicesAE,test{2});
%% Calculate the volume
stats.VolumeAE = calculateVolume(stats.FeaturesAE);
%% Calculate the significancy for volume
stats.sVolumeAE = calculateSvolume(stats.VolumeAE,test{2});
%% Calculate the change in volume
stats.cVolumeAE = calculateCvolume(stats.FeaturesAE);
%% Calculate significancy forthe change in volume
stats.sCvolumeAE = calculateSCvolume(stats.cVolumeAE,test{2});

%---- Compare between before and after implementing QCS
%% Calculate the significancy for ICC between before and after
stats.sICCBA = calculateSiccBA(stats.ICC,stats.ICCAE);
%% Calculate the significancy for ICC between before and after (qith random 1000)
stats.sICCBAR = calculateSiccBAR(stats.ICCR,stats.ICCAE);

%---- Before including post-processing step
BPP = load('features-beforePP.mat');
esbppPath = fullfile(masterPath,'Data','Segmentation','Automated');
%% Calculate ICC
stats.ICCBPP = calculateICC(BPP.features);
%% Calculate the significancy for ICC
stats.sICCBPP = calculateSicc(stats.ICCBPP);
%% Calculate CV
stats.CVBPP = calculateCV(BPP.features);
%% Calculate the significancy for CV
stats.sCVBPP = calculateScv(stats.CVBPP,test{1});
%% Calculate the perfromance scores
stats.pScoresBPP = qualityControl(masterPath,esbppPath,factors,0);
%% Calculate the significancy for performance metrics
stats.sScoresBPP = calculateSscores(stats.pScoresBPP.pScores,test{1});
%% Calculate number of included slices
stats.nSlicesBPP = calculateNslices(masterPath,esbppPath,[]);
%% Calculate the significancy for number of included slices
stats.sSlicesBPP = calculateSslices(stats.nSlicesBPP,test{1});
%% Calculate the volume
stats.VolumeBPP = calculateVolume(BPP.features);
%% Calculate the significancy for volume
stats.sVolumeBPP = calculateSvolume(stats.VolumeBPP,test{1});
%% Calculate the change in volume
stats.cVolumeBPP = calculateCvolume(BPP.features);
%% Calculate significancy forthe change in volume
stats.sCvolumeBPP = calculateSCvolume(stats.cVolumeBPP,test{1});

%---- Compare between before and after including post-processing
%% Calculate the significancy for ICC between before and after post-process step
stats.sICCBABPP = calculateSiccBA(stats.ICC,stats.ICCBPP);

%% Save
save('stats.mat','stats')
end

%% getProtocol
% This function is to get the scanning protocol parameters for Scan 1 and
% Scan 2 cases
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%
% output:
%   1. protocol: The scanning protocol parameters. (structure)
%
function protocol = getProtocol(masterPath)
% DICOM folders path
dicomP = '.../dataset';
% Original Cases path
oPath = fullfile(masterPath,'Data','Original');
% Loop over Scans
oD = dir(oPath);
oD = oD(~ismember({oD.name},{'.','..'}));
for kk = 1:numel(oD)
    % Make table to fill later
    protocol.(oD(kk).name) = table;
    % Loop over cases
    cD = dir(fullfile(oD(kk).folder,oD(kk).name,'*.mhd'));
    for cc = 1:numel(cD)
        % Name
        protocol.(oD(kk).name).Name(cc,:) = {cD(cc).name(1:18)};
        % Get into the case DICOM folder
        dicomD = dir(fullfile(dicomP,cD(cc).name(1:8),...
            cD(cc).name(1:18)));
        % Find the correct T2W secquance
        % Patterns we want
        pattern = ["T2W_tra","T2W tra","T2W TRA","T2W_TRA",...
            "t2_tse_tra","t2 tse tra","t2 tse _ra","t2_tse tra"];
        % Secquances with the pattern
        sIdx = contains({dicomD.name},pattern,'IgnoreCase',true);
        imagefile = dicomD(sIdx);
        % In case more than one match, order them
        infoTemp = zeros(1,numel(imagefile));
        for mm = 1:numel(imagefile)
            sD = dir(fullfile(imagefile(mm).folder,imagefile(mm).name));
            sD = sD(~ismember({sD.name},{'.','..'}));
            inf = dicominfo(fullfile(sD(1).folder,sD(1).name));
            infoTemp(mm) = inf.SeriesNumber;
        end
        [~,maxSIdx] = max(infoTemp);
        [~,minSIdx] = min(infoTemp);
        % Read the newest sequance except if it was an exception
        exceptions = {'.....'};
        if ismember(cD(cc).name(1:18),exceptions)
            seq = fullfile(imagefile(minSIdx).folder,imagefile(minSIdx).name);% espcial cases
        else
            seq = fullfile(imagefile(maxSIdx).folder,imagefile(maxSIdx).name);
        end
        % Get the dicom tags
        iD = dir(seq);
        iD = iD(~ismember({iD.name},{'.','..'}));
        info = dicominfo(fullfile(iD(1).folder,iD(1).name));
        
        protocol.(oD(kk).name).PatientAge(cc,:) = str2double(info.PatientAge(1:end-1));
        protocol.(oD(kk).name).ScanDate(cc,:) =...
            datestr(datenum(num2str(str2double(info.SeriesDate),'%d'),'yyyymmdd'),'yyyy/mm/dd');
        protocol.(oD(kk).name).MagneticFieldStrength(cc,:) = info.MagneticFieldStrength;
        protocol.(oD(kk).name).Manufacturer(cc,:) = {info.Manufacturer};
        protocol.(oD(kk).name).ManufacturerModelName(cc,:) = {info.ManufacturerModelName};
        protocol.(oD(kk).name).ScanningSequence(cc,:) = {info.ScanningSequence};
        protocol.(oD(kk).name).RepetitionTime(cc,:) = info.RepetitionTime;
        protocol.(oD(kk).name).EchoTime(cc,:) = info.EchoTime;
        protocol.(oD(kk).name).FlipAngle(cc,:) = info.FlipAngle;
        protocol.(oD(kk).name).NumberOfAverages(cc,:) = info.NumberOfAverages;
        protocol.(oD(kk).name).SlicesNumber(cc,:) = numel(iD)-2;
        protocol.(oD(kk).name).Width(cc,:) = info.Width;
        protocol.(oD(kk).name).Height(cc,:) = info.Height;
        protocol.(oD(kk).name).SliceThickness(cc,:) = info.SliceThickness;
        protocol.(oD(kk).name).PixelSpacing(cc,:) = info.PixelSpacing(1);
    end
end
end

%% getProtocolTr
% This function is to get the scanning protocol parameters for segmentation
% training set.
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%
% output:
%   1. protocol: The scanning protocol parameters. (table)
%
function protocol = getProtocolTr(masterPath)
% DICOM folders path
dicomP = '.../dataset';
% Original Cases path
oPath = fullfile(masterPath,'Data','Train','Original');

% Make table to fill later
protocol = table;
% Loop over cases
cD = dir(fullfile(oPath,'*.mhd'));
for cc = 1:numel(cD)
    % Name
    protocol.Name(cc,:) = {cD(cc).name(1:18)};
    % Get into the case DICOM folder
    dicomD = dir(fullfile(dicomP,cD(cc).name(1:8),...
        cD(cc).name(1:18)));
    % Find the correct T2W secquance
    % Patterns we want
    pattern = ["T2W_tra","T2W tra","T2W TRA","T2W_TRA",...
        "t2_tse_tra","t2 tse tra","t2 tse _ra","t2_tse tra"];
    % Secquances with the pattern
    sIdx = contains({dicomD.name},pattern,'IgnoreCase',true);
    imagefile = dicomD(sIdx);
    % In case more than one match, order them
    infoTemp = zeros(1,numel(imagefile));
    for mm = 1:numel(imagefile)
        sD = dir(fullfile(imagefile(mm).folder,imagefile(mm).name));
        sD = sD(~ismember({sD.name},{'.','..'}));
        inf = dicominfo(fullfile(sD(1).folder,sD(1).name));
        infoTemp(mm) = inf.SeriesNumber;
    end
    [~,maxSIdx] = max(infoTemp);
    [~,minSIdx] = min(infoTemp);
    % Read the newest sequance except if it was an exception
    exceptions = {'.....'};
    if ismember(cD(cc).name(1:18),exceptions)
        seq = fullfile(imagefile(minSIdx).folder,imagefile(minSIdx).name);% espcial cases
    else
        seq = fullfile(imagefile(maxSIdx).folder,imagefile(maxSIdx).name);
    end
    % Get the dicom tags
    iD = dir(seq);
    iD = iD(~ismember({iD.name},{'.','..'}));
    info = dicominfo(fullfile(iD(1).folder,iD(1).name));
    
    protocol.PatientAge(cc,:) = str2double(info.PatientAge(1:end-1));
    protocol.ScanDate(cc,:) =...
        datestr(datenum(num2str(str2double(info.SeriesDate),'%d'),'yyyymmdd'),'yyyy/mm/dd');
    protocol.MagneticFieldStrength(cc,:) = info.MagneticFieldStrength;
    protocol.Manufacturer(cc,:) = {info.Manufacturer};
    protocol.ManufacturerModelName(cc,:) = {info.ManufacturerModelName};
    protocol.ScanningSequence(cc,:) = {info.ScanningSequence};
    protocol.RepetitionTime(cc,:) = info.RepetitionTime;
    protocol.EchoTime(cc,:) = info.EchoTime;
    protocol.FlipAngle(cc,:) = info.FlipAngle;
    protocol.NumberOfAverages(cc,:) = info.NumberOfAverages;
    protocol.SlicesNumber(cc,:) = numel(iD)-2;
    protocol.Width(cc,:) = info.Width;
    protocol.Height(cc,:) = info.Height;
    protocol.SliceThickness(cc,:) = info.SliceThickness;
    protocol.PixelSpacing(cc,:) = info.PixelSpacing(1);
end
end

%% calculateICC
% This function is to calculate the intraclass correlation (ICC)
% between the features those extacted from Scan 1 and Scan 2
%
% Input:
%   1. features: The extracted shape features. (structure)
%
% output:
%   1. ICC: The calculated ICC results. (structure)
%
function ICC = calculateICC(features)
% Get Netwroks names
nets = fieldnames(features);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions names
    regs = fieldnames(features.(nets{ii}));
    % Make tables to fill
    ICC.(nets{ii}).ICC = table;
    ICC.(nets{ii}).ICCpVal = table;
    ICC.(nets{ii}).ICCCI = table;
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scans names
        scans = fieldnames(features.(nets{ii}).(regs{jj}));
        % Get features names
        fNames = features.V_Net_3D.PZ.Scan1.Properties.VariableNames;
        % Loop over features to calculate ICC
        for ff = 1:numel(fNames)
            % Dummy matrix
            temp = [features.(nets{ii}).(regs{jj}).(scans{1}).(fNames{ff}),...
                features.(nets{ii}).(regs{jj}).(scans{2}).(fNames{ff})];
            % Set alpha
            alpha = 0.05;
            % Temp ICC
            tICC = f_ICC(temp,alpha);
            % Assign to table
            ICC.(nets{ii}).ICC.(regs{jj})(ff,1) = tICC{1, 2}.est;
            ICC.(nets{ii}).ICCpVal.(regs{jj})(ff,1) = tICC{1, 2}.pVal;
            ICC.(nets{ii}).ICCCI.(regs{jj}){ff,1} = tICC{1, 2}.confInterval;
            % clear
            clear temp tICC
        end
    end
    % Set features names as raws
    rnames = fNames';
    ICC.(nets{ii}).ICC.Properties.RowNames = regexprep(rnames, 'shape_', '');
    ICC.(nets{ii}).ICCpVal.Properties.RowNames = regexprep(rnames, 'shape_', '');
    ICC.(nets{ii}).ICCCI.Properties.RowNames = regexprep(rnames, 'shape_', '');
end
end

%% calculateSicc
% This function is to calculate the significance of changes of the
% calculated ICC
%
% Input:
%   1. ICC: The calculated ICC results. (structure)
%
% output:
%   1. sICC: The significance of ICC changes. (structure)
%
function sICC = calculateSicc(ICC)
% --Test between Networks for All features
% Get Netwroks names
nets = fieldnames(ICC);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Make table to fill later
    sICC.NetsAll = table;
    % Get regions names
    regs = ICC.(nets{ii}).ICC.Properties.VariableNames;
    % Loop over regions
    for jj = 1:numel(regs)
        % Do the statistical test
        for ll = 1:numel(nets)
            [pvNA.(regs{jj})(ll,ii),~] =...
                signrank(ICC.(nets{ll}).ICC.(regs{jj})(:),...
                ICC.(nets{ii}).ICC.(regs{jj})(:));
        end
    end
end

% Get raw names
netsCN = cell(numel(netsN),numel(netsN));
for ii = 1:numel(netsN)
    for jj = ii:numel(netsN)-1
        netsCN{jj,ii} = [netsN{ii} ' Vs ' netsN{jj+1}];
    end
end
netsCN = netsCN(~cellfun('isempty',netsCN));
netsCN = netsCN(:);

% Loop over regions to correct
for jj = 1:numel(regs)
    % Correct
    inP = nonzeros(tril(pvNA.(regs{jj}),-1));
    [~, ~, ~, sICC.NetsAll.(regs{jj})] =...
        fdr_bh(inP,0.05,'pdep','no');
end
% Add raw names
sICC.NetsAll.Properties.RowNames = netsCN;

% --Test between Regions
% Loop over Netwroks
for ii = 1:numel(nets)
    % Make table to fill later
    sICC.Regions = table;
    % Loop over regions
    for jj = 1:numel(regs)
        % Do the statistical test
        for ll = 1:numel(regs)
            [pvR.(nets{ii})(ll,jj),~] =...
                signrank(ICC.(nets{ii}).ICC.(regs{jj})(:),...
                ICC.(nets{ii}).ICC.(regs{ll})(:));
        end
    end
end

% Get raw names
reCN = cell(numel(regs),numel(regs));
for ii = 1:numel(regs)
    for jj = ii:numel(regs)-1
        reCN{jj,ii} = [regs{ii} ' Vs ' regs{jj+1}];
    end
end
reCN = reCN(~cellfun('isempty',reCN));
reCN = reCN(:);

% Loop over Netwroks
for ii = 1:numel(nets)
    % Correct
    inP = nonzeros(tril(pvR.(nets{ii}),-1));
    [~, ~, ~, sICC.Regions.(nets{ii})] =...
        fdr_bh(inP,0.05,'pdep','no');
end
% Add raw names
sICC.Regions.Properties.RowNames = reCN;

% --Test between Networks per feature
% Loop over Netwroks
for ii = 1:numel(nets)
    % Make table to fill later
    sICC.Nets.(nets{ii}) = table;
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over features
        fNames = ICC.(nets{ii}).ICCCI.Properties.RowNames;
        for ff = 1:numel(fNames)
            % Do the statistical test
            sICC.Nets.(nets{ii}).(regs{jj})(ff) =...
                isempty(range_intersection(ICC.(nets{1}).ICCCI.(regs{jj}){ff},...
                ICC.(nets{ii}).ICCCI.(regs{jj}){ff}));
        end
    end
end
end

%% calculateCV
% This function is to calculate the coefficent of variation (CV)
% between the features those extacted from Scan 1 and Scan 2
%
% Input:
%   1. features: The extracted shape features. (structure)
%
% output:
%   1. CV: The calculated CV results. (structure)
%
function CV = calculateCV(features)
% Get Netwroks names
nets = fieldnames(features);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions names
    regs = fieldnames(features.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scans names
        scans = fieldnames(features.(nets{ii}).(regs{jj}));
        % Get features names
        fNames = features.V_Net_3D.PZ.Scan1.Properties.VariableNames;
        % Loop over features to calculate ICC
        for ff = 1:numel(fNames)
            % Loop over cases
            for ll = 1:size(features.(nets{ii}).(regs{jj}).(scans{1}).(fNames{ff}),1)
                % values
                v1 = features.(nets{ii}).(regs{jj}).(scans{1}).(fNames{ff})(ll);
                v2 = features.(nets{ii}).(regs{jj}).(scans{2}).(fNames{ff})(ll);
                CV.(nets{ii}).(regs{jj})(ll,ff) = (std([v1,v2]))/(mean([v1,v2]));
            end
        end
        % Convert to table
        CV.(nets{ii}).(regs{jj}) = array2table(CV.(nets{ii}).(regs{jj}));
        % Set features names as columns
        vnames = fNames';
        CV.(nets{ii}).(regs{jj}).Properties.VariableNames = regexprep(vnames, 'shape_', '');
        % Set cases as raws
        rname = features.(nets{ii}).(regs{jj}).(scans{1}).Properties.RowNames;
        rname = cellfun(@(x) x(1:8), rname, 'un', 0);
        CV.(nets{ii}).(regs{jj}).Properties.RowNames = rname;
    end
end
end

%% calculateScv
% This function is to calculate the significance of changes of the
% calculated CV
%
% Input:
%   1. CV: The calculated CV results. (structure)
%   2. test: the statistical test. (string)
%
% output:
%   1. sCV: The significance of CV changes. (structure)
%
function sCV = calculateScv(CV,test)
% --Test between Networks
% Get Netwroks names
nets = fieldnames(CV);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(CV.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Make table to fill later
        sCV.Nets.(regs{jj}) = table;
        % Get features names
        fnames = CV.(nets{ii}).(regs{jj}).Properties.VariableNames;
        % Do the statistical test
        for ff = 1:numel(fnames)
            for ll = 1:numel(nets)
                if contains(test,'signrank')
                    [pvN.(regs{jj}).(fnames{ff})(ll,ii),~] =...
                        signrank(table2array(CV.(nets{ll}).(regs{jj})(:,ff)),...
                        table2array(CV.(nets{ii}).(regs{jj})(:,ff)));
                elseif contains(test,'ranksum')
                    [pvN.(regs{jj}).(fnames{ff})(ll,ii),~] =...
                        ranksum(table2array(CV.(nets{ll}).(regs{jj})(:,ff)),...
                        table2array(CV.(nets{ii}).(regs{jj})(:,ff)));
                end
            end
        end
    end
end

% Get raw names
netsCN = cell(numel(netsN),numel(netsN));
for ii = 1:numel(netsN)
    for jj = ii:numel(netsN)-1
        netsCN{jj,ii} = [netsN{ii} ' Vs ' netsN{jj+1}];
    end
end
netsCN = netsCN(~cellfun('isempty',netsCN));
netsCN = netsCN(:);

% Loop over regions to correct
for jj = 1:numel(regs)
    % Loop over features
    for ff = 1:numel(fnames)
        % Correct
        inP = nonzeros(tril(pvN.(regs{jj}).(fnames{ff}),-1));
        [~, ~, ~, sCV.Nets.(regs{jj}).(fnames{ff})] =...
            fdr_bh(inP,0.05,'pdep','no');
    end
    % Add raw names
    sCV.Nets.(regs{jj}).Properties.RowNames = netsCN;
end

% --Test between features
% Loop over Netwroks
for ii = 1:numel(nets)
    % Make table to fill later
    sCV.Features.(nets{ii}) = table;
    % Loop over regions
    for jj = 1:numel(regs)
        % Do the statistical test
        for ff = 1:numel(fnames)
            for ll = 1:numel(fnames)
                if contains(test,'signrank')
                    [pvF.(nets{ii}).(regs{jj})(ll,ff),~] =...
                        signrank(table2array(CV.(nets{ii}).(regs{jj})(:,ff)),...
                        table2array(CV.(nets{ii}).(regs{jj})(:,ll)));
                elseif contains(test,'ranksum')
                    [pvF.(nets{ii}).(regs{jj})(ll,ff),~] =...
                        ranksum(table2array(CV.(nets{ii}).(regs{jj})(:,ff)),...
                        table2array(CV.(nets{ii}).(regs{jj})(:,ll)));
                end
            end
        end
    end
end

% Get raw names
feCN = cell(numel(fnames),numel(fnames));
for ii = 1:numel(fnames)
    for jj = ii:numel(fnames)-1
        feCN{jj,ii} = [fnames{ii} ' Vs ' fnames{jj+1}];
    end
end
feCN = feCN(~cellfun('isempty',feCN));
feCN = feCN(:);

% Loop over Netwroks
for ii = 1:numel(nets)
    % Loop over regions to correct
    for jj = 1:numel(regs)
        % Correct
        for ff = 1:numel(fnames)
            inP = nonzeros(tril(pvF.(nets{ii}).(regs{jj}),-1));
            [~, ~, ~, sCV.Features.(nets{ii}).(regs{jj})] =...
                fdr_bh(inP,0.05,'pdep','no');
        end
        % Add raw names
        sCV.Features.(nets{ii}).Properties.RowNames = feCN;
    end
end

% --Test between regions
% Loop over Netwroks
for ii = 1:numel(nets)
    % Make table to fill later
    sCV.Regions.(nets{ii}) = table;
    % Loop over regions
    for jj = 1:numel(regs)
        % Do the statistical test
        for ff = 1:numel(fnames)
            for ll = 1:numel(regs)
                if contains(test,'signrank')
                    [pvR.(nets{ii}).(fnames{ff})(ll,jj),~] =...
                        signrank(table2array(CV.(nets{ii}).(regs{jj})(:,ff)),...
                        table2array(CV.(nets{ii}).(regs{ll})(:,ff)));
                elseif contains(test,'ranksum')
                    [pvR.(nets{ii}).(fnames{ff})(ll,jj),~] =...
                        ranksum(table2array(CV.(nets{ii}).(regs{jj})(:,ff)),...
                        table2array(CV.(nets{ii}).(regs{ll})(:,ff)));
                end
            end
        end
    end
end

% Get raw names
reCN = cell(numel(regs),numel(regs));
for ii = 1:numel(regs)
    for jj = ii:numel(regs)-1
        reCN{jj,ii} = [regs{ii} ' Vs ' regs{jj+1}];
    end
end
reCN = reCN(~cellfun('isempty',reCN));
reCN = reCN(:);

% Loop over Netwroks
for ii = 1:numel(nets)
    % Loop over regions to correct
    for jj = 1:numel(regs)
        % Correct
        for ff = 1:numel(fnames)
            inP = nonzeros(tril(pvR.(nets{ii}).(fnames{ff}),-1));
            [~, ~, ~, sCV.Regions.(nets{ii}).(fnames{ff})] =...
                fdr_bh(inP,0.05,'pdep','no');
        end
        % Add raw names
        sCV.Regions.(nets{ii}).Properties.RowNames = reCN;
    end
end
end

%% calculateSscores
% This function is to calculate the significance of changes of the
% calculated performance scores
%
% Input:
%   1. scores: The calculated performance scores for the cases. (structure)
%   2. test: the statistical test. (string)
%
% output:
%   1. sCV: The significance of the changes in the perfrormance scores. (structure)
%
function sScores = calculateSscores(scores,test)
% --Test between Scan 1 and Scan 2
% Get Netwroks names
nets = fieldnames(scores);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(scores.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scan names
        scans = fieldnames(scores.(nets{ii}).(regs{jj}));
        % Get the metrics fields
        mets = fieldnames(scores.(nets{ii}).(regs{jj}).(scans{1}));
        mets = mets(~contains(mets,'names'));
        % Loop over the metrics
        for mm = 1:numel(mets)
            % Check if it is structure
            if isstruct(scores.(nets{ii}).(regs{jj}).(scans{1}).(mets{mm}))
                % Get metrics from that field
                mf = fieldnames(scores.(nets{ii}).(regs{jj}).(scans{1}).(mets{mm}));
                for ff = 1:numel(mf)
                    % Make table to fill
                    sScores.Scans.(mets{mm}).(mf{ff}) = table;
                    % Do the statistical test
                    if contains(test,'signrank')
                        [pv.(mets{mm}).(mf{ff}).(regs{jj})(ii,:),~] = ...
                            signrank(scores.(nets{ii}).(regs{jj}).(scans{1}).(mets{mm}).(mf{ff})(:),...
                            scores.(nets{ii}).(regs{jj}).(scans{2}).(mets{mm}).(mf{ff})(:));
                    elseif contains(test,'ranksum')
                        [pv.(mets{mm}).(mf{ff}).(regs{jj})(ii,:),~] = ...
                            ranksum(scores.(nets{ii}).(regs{jj}).(scans{1}).(mets{mm}).(mf{ff})(:),...
                            scores.(nets{ii}).(regs{jj}).(scans{2}).(mets{mm}).(mf{ff})(:));
                    end
                end
            else
                % Make table to fill
                sScores.Scans.(mets{mm}) = table;
                % Do the statistical test
                [pv.(mets{mm}).(regs{jj})(ii,:),~] =...
                    signrank(scores.(nets{ii}).(regs{jj}).(scans{1}).(mets{mm})(:),...
                    scores.(nets{ii}).(regs{jj}).(scans{2}).(mets{mm})(:));
            end
        end
    end
end

% Loop over regions
for jj = 1:numel(regs)
    % Loop over the metrics
    for mm = 1:numel(mets)
        % Check if it is structure
        if isstruct(scores.(nets{ii}).(regs{jj}).(scans{1}).(mets{mm}))
            % Loop over metrics from that field
            for ff = 1:numel(mf)
                % Correct
                [~, ~, ~,  sScores.Scans.(mets{mm}).(mf{ff}).(regs{jj})] =...
                    fdr_bh(pv.(mets{mm}).(mf{ff}).(regs{jj}),0.05,'pdep','no');
                % Add raw names
                sScores.Scans.(mets{mm}).(mf{ff}).Properties.RowNames = netsN;
            end
        else
            % Correct
            [~, ~, ~, sScores.Scans.(mets{mm}).(regs{jj})] =...
                fdr_bh(pv.(mets{mm}).(regs{jj}),0.05,'pdep','no');
            % Add raw names
            sScores.Scans.(mets{mm}).Properties.RowNames = netsN;
        end
    end
end

% --Test between netwroks
% Loop over Netwroks
for ii = 1:numel(nets)
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over scans
        for kk = 1:numel(scans)
            % Loop over the metrics
            for mm = 1:numel(mets)
                % Check if it is structure
                if isstruct(scores.(nets{ii}).(regs{jj}).(scans{1}).(mets{mm}))
                    % Loop over the metrics from that field
                    for ff = 1:numel(mf)
                        % Make table to fill
                        sScores.(['Nets' (scans{kk})]).(mets{mm}).(mf{ff}) = table;
                        % Do the statistical test
                        for ll = 1:numel(nets)
                            if contains(test,'signrank')
                                [pvN.(mets{mm}).(mf{ff}).(regs{jj}).(scans{kk})(ll,ii),~] =....
                                    signrank(scores.(nets{ll}).(regs{jj}).(scans{kk}).(mets{mm}).(mf{ff})(:),...
                                    scores.(nets{ii}).(regs{jj}).(scans{kk}).(mets{mm}).(mf{ff})(:));
                            elseif contains(test,'ranksum')
                                [pvN.(mets{mm}).(mf{ff}).(regs{jj}).(scans{kk})(ll,ii),~] =....
                                    ranksum(scores.(nets{ll}).(regs{jj}).(scans{kk}).(mets{mm}).(mf{ff})(:),...
                                    scores.(nets{ii}).(regs{jj}).(scans{kk}).(mets{mm}).(mf{ff})(:));
                            end
                        end
                    end
                else
                    % Make table to fill
                    sScores.(['Nets' (scans{kk})]).(mets{mm}) = table;
                    % Do the statistical test
                    for ll = 1:numel(nets)
                        if contains(test,'signrank')
                            [pvN.(mets{mm}).(regs{jj}).(scans{kk})(ll,ii),~] =....
                                signrank(scores.(nets{ll}).(regs{jj}).(scans{kk}).(mets{mm})(:),...
                                scores.(nets{ii}).(regs{jj}).(scans{kk}).(mets{mm})(:));
                        elseif contains(test,'ranksum')
                            [pvN.(mets{mm}).(regs{jj}).(scans{kk})(ll,ii),~] =....
                                ranksum(scores.(nets{ll}).(regs{jj}).(scans{kk}).(mets{mm})(:),...
                                scores.(nets{ii}).(regs{jj}).(scans{kk}).(mets{mm})(:));
                        end
                    end
                end
            end
        end
    end
end

% Get raw names
netsCN = cell(numel(netsN),numel(netsN));
for ii = 1:numel(netsN)
    for jj = ii:numel(netsN)-1
        netsCN{jj,ii} = [netsN{ii} ' Vs ' netsN{jj+1}];
    end
end
netsCN = netsCN(~cellfun('isempty',netsCN));
netsCN = netsCN(:);

% Loop over regions
for jj = 1:numel(regs)
    % Loop over the metrics
    for mm = 1:numel(mets)
        % Check if it is structure
        if isstruct(scores.(nets{ii}).(regs{jj}).(scans{1}).(mets{mm}))
            % Loop over metrics from that field
            for ff = 1:numel(mf)
                % Loop over scans
                for kk = 1:numel(scans)
                    % Correct
                    inP = nonzeros(tril(pvN.(mets{mm}).(mf{ff}).(regs{jj}).(scans{kk}),-1));
                    [~, ~, ~, sScores.(['Nets' (scans{kk})]).(mets{mm}).(mf{ff}).(regs{jj})] =...
                        fdr_bh(inP,0.05,'pdep','no');
                    % Add raw names
                    sScores.(['Nets' (scans{kk})]).(mets{mm}).(mf{ff}).Properties.RowNames = netsCN;
                end
            end
        else
            % Loop over scans
            for kk = 1:numel(scans)
                % Correct
                inP = nonzeros(tril(pvN.(mets{mm}).(regs{jj}).(scans{kk}),-1));
                [~, ~, ~,sScores.(['Nets' (scans{kk})]).(mets{mm}).(regs{jj})] =...
                    fdr_bh(inP,0.05,'pdep','no');
                % Add raw names
                sScores.(['Nets' (scans{kk})]).(mets{mm}).Properties.RowNames = netsCN;
            end
        end
    end
end
end

%% calculateNslices
% This function is to calculate the number of the included slices in
% calculations for the cases in Scan1 and Scan2
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%   2. esPath: The path to the automatically segmented cases post processing. (string)
%   3. eIdx: The excluded cases Idxs. (structure)
%
% output:
%   1. nSlices: The numbers of the included slices. (structure)
%
function nSlices = calculateNslices(masterPath,esPath,eIdx)
% Paths
msPath = fullfile(masterPath,'Data','Segmentation','Manual','Final');
% Loop over networks
taD = dir(esPath);
taD = taD(~ismember({taD.name},{'.','..'}));
aD(1).name = 'Manual'; % add manual
aD(1).folder = msPath; % add manual path
aD(1).date = taD(1).date;
aD(1).bytes = taD(1).bytes;
aD(1).isdir = taD(1).isdir;
aD(1).datenum = taD(1).datenum;
aD(2:numel(taD)+1) = taD(1:end);
for ii = 1:numel(aD)
    % Get regions names
    if strcmp(aD(ii).name,'Manual')
        nD = dir(fullfile(esPath,aD(ii+1).name));
        nD = nD(~ismember({nD.name},{'.','..'}));
    else
        nD = dir(fullfile(esPath,aD(ii).name));
        nD = nD(~ismember({nD.name},{'.','..'}));
    end
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
            
            if strcmp(aD(ii).name,'Manual')
                sD = dir(fullfile(aD(ii).folder,rD(kk).name,nD(jj).name,'*.mhd'));
            else
                sD = dir(fullfile(rD(kk).folder,rD(kk).name,'*.mhd'));
                if ~isempty(eIdx)
                    sD(eIdx.(netName)) = [];
                end
            end
            
            for ll = 1:numel(sD)
                FilePath = fullfile(sD(ll).folder,sD(ll).name);
                nSlices.(netName).(nD(jj).name).(rD(kk).name)(ll,:) =...
                    getNslices(FilePath);
            end
        end
    end
end
end

%% getNslices
% This function is to calculate the number of the included slices in a mask
%
% Input:
%   1. FilePath: The path to the case. (string)
%
% output:
%   1. nSlices: The numbers of the included slices. (structure)
%
function nSlices = getNslices(FilePath)
% Read the image in matlab
[StrDatax, ~, ~] = elxMetaIOFileToStrDatax(FilePath, 0);
% Loop over all slices and find if it contains mask
sIdx = zeros(size(StrDatax.Data,3),1);
for ii = 1:size(StrDatax.Data,3)
    sIdx(ii,:) = sum(sum(StrDatax.Data(:,:,ii)))>0;
end
nSlices = sum(sIdx);
end

%% calculateSslices
% This function is to calculate the significance of change between
% Scan 1 and Scan 2 volumes
%
% Input:
%   1. nSlices: The numbers of the included slices. (structure)
%   2. test: the statistical test. (string)
%
% output:
%   1. sSlices: The significance of slices numbers change. (structure)
%
function sSlices = calculateSslices(nSlices,test)
% --Test between Scan 1 and Scan 2
sSlices.Scans = table;
% Get Netwroks names
nets = fieldnames(nSlices);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(nSlices.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scan names
        scans = fieldnames(nSlices.(nets{ii}).(regs{jj}));
        % Do the statistical test
        if contains(test,'signrank')
            [pv.(regs{jj})(ii,:),~] =...
                signrank(nSlices.(nets{ii}).(regs{jj}).(scans{1})(:),...
                nSlices.(nets{ii}).(regs{jj}).(scans{2})(:));
        elseif contains(test,'ranksum')
            [pv.(regs{jj})(ii,:),~] =...
                ranksum(nSlices.(nets{ii}).(regs{jj}).(scans{1})(:),...
                nSlices.(nets{ii}).(regs{jj}).(scans{2})(:));
        end
    end
end
% Loop over regions
for jj = 1:numel(regs)
    % Correct
    [~, ~, ~, sSlices.Scans.(regs{jj})] = fdr_bh(pv.(regs{jj}),0.05,'pdep','no');
end
% Add raw names
sSlices.Scans.Properties.RowNames = netsN;

% --Test between netwroks
% Loop over Netwroks
for ii = 1:numel(nets)
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over scans
        for kk = 1:numel(scans)
            sSlices.(['Nets' (scans{kk})]) = table;
            % Do the statistical test
            for ll = 1:numel(nets)
                if contains(test,'signrank')
                    [pvN.(regs{jj}).(scans{kk})(ll,ii),~] =...
                        signrank(nSlices.(nets{ll}).(regs{jj}).(scans{kk})(:),...
                        nSlices.(nets{ii}).(regs{jj}).(scans{kk})(:));
                elseif contains(test,'ranksum')
                    [pvN.(regs{jj}).(scans{kk})(ll,ii),~] =...
                        ranksum(nSlices.(nets{ll}).(regs{jj}).(scans{kk})(:),...
                        nSlices.(nets{ii}).(regs{jj}).(scans{kk})(:));
                end
            end
        end
    end
end

% Loop over regions
for jj = 1:numel(regs)
    % Loop over scans
    for kk = 1:numel(scans)
        % Correct
        inP = nonzeros(tril(pvN.(regs{jj}).(scans{kk}),-1));
        [~, ~, ~, sSlices.(['Nets' (scans{kk})]).(regs{jj})] =...
            fdr_bh(inP,0.05,'pdep','no');
    end
end
% Add raw names
netsCN = cell(numel(netsN),numel(netsN));
for ii = 1:numel(netsN)
    for jj = ii:numel(netsN)-1
        netsCN{jj,ii} = [netsN{ii} ' Vs ' netsN{jj+1}];
    end
end
netsCN = netsCN(~cellfun('isempty',netsCN));
netsCN = netsCN(:);

sSlices.NetsScan1.Properties.RowNames = netsCN;
sSlices.NetsScan2.Properties.RowNames = netsCN;
end

%% calculateVolume
% This function is to calculate the Scan 1 and Scan 2 volumes
%
% Input:
%   1. features: The extracted shape features. (structure)
%
% output:
%   1. Volume: The volumes. (structure)
%
function Volume = calculateVolume(features)
% Get Netwroks names
nets = fieldnames(features);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions names
    regs = fieldnames(features.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over scans
        scans = fieldnames(features.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            % Loop over cases
            for ll = 1:size(features.(nets{ii}).(regs{jj}).(scans{kk}),1)
                Volume.(nets{ii}).(regs{jj}).(scans{kk})(ll,:) =...
                    table2array(features.(nets{ii}).(regs{jj}).(scans{kk})(ll,14));
            end
        end
    end
end
end

%% calculateSvolume
% This function is to calculate the significance of change between
% Scan 1 and Scan 2 volumes
%
% Input:
%   1. Volume: The volumes. (structure)
%   2. test: the statistical test. (string)
%
% output:
%   1. sVolume: The significance of volumes change. (structure)
%
function sVolume = calculateSvolume(Volume,test)
% --Test between Scan 1 and Scan 2
sVolume.Scans = table;
% Get Netwroks names
nets = fieldnames(Volume);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(Volume.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scan names
        scans = fieldnames(Volume.(nets{ii}).(regs{jj}));
        % Do the statistical test
        if contains(test,'signrank')
            [pv.(regs{jj})(ii,:),~] =...
                signrank(Volume.(nets{ii}).(regs{jj}).(scans{1})(:),...
                Volume.(nets{ii}).(regs{jj}).(scans{2})(:));
        elseif contains(test,'ranksum')
            [pv.(regs{jj})(ii,:),~] =...
                ranksum(Volume.(nets{ii}).(regs{jj}).(scans{1})(:),...
                Volume.(nets{ii}).(regs{jj}).(scans{2})(:));
        end
    end
end

% Loop over regions
for jj = 1:numel(regs)
    % Correct
    [~, ~, ~, sVolume.Scans.(regs{jj})] =...
        fdr_bh(pv.(regs{jj}),0.05,'pdep','no');
end
% Add raw names
sVolume.Scans.Properties.RowNames = netsN;

% --Test between netwroks
% Loop over Netwroks
for ii = 1:numel(nets)
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over scans
        for kk = 1:numel(scans)
            sVolume.(['Nets' (scans{kk})]) = table;
            % Do the statistical test
            for ll = 1:numel(nets)
                if contains(test,'signrank')
                    [pvN.(regs{jj}).(scans{kk})(ll,ii),~] =...
                        signrank(Volume.(nets{ll}).(regs{jj}).(scans{kk})(:),...
                        Volume.(nets{ii}).(regs{jj}).(scans{kk})(:));
                elseif contains(test,'ranksum')
                    [pvN.(regs{jj}).(scans{kk})(ll,ii),~] =...
                        ranksum(Volume.(nets{ll}).(regs{jj}).(scans{kk})(:),...
                        Volume.(nets{ii}).(regs{jj}).(scans{kk})(:));
                end
            end
        end
    end
end

% Loop over regions
for jj = 1:numel(regs)
    % Loop over scans
    for kk = 1:numel(scans)
        % Correct
        inP = nonzeros(tril(pvN.(regs{jj}).(scans{kk}),-1));
        [~, ~, ~, sVolume.(['Nets' (scans{kk})]).(regs{jj})] =...
            fdr_bh(inP,0.05,'pdep','no');
    end
end

% Add raw names
netsCN = cell(numel(netsN),numel(netsN));
for ii = 1:numel(netsN)
    for jj = ii:numel(netsN)-1
        netsCN{jj,ii} = [netsN{ii} ' Vs ' netsN{jj+1}];
    end
end
netsCN = netsCN(~cellfun('isempty',netsCN));
netsCN = netsCN(:);

sVolume.NetsScan1.Properties.RowNames = netsCN;
sVolume.NetsScan2.Properties.RowNames = netsCN;
end

%% calculateCvolume
% This function is to calculate the percentage of change between
% Scan 1 and Scan 2 volumes
%
% Input:
%   1. features: The extracted shape features. (structure)
%
% output:
%   1. cVolume: The change in volumes. (structure)
%
function cVolume = calculateCvolume(features)
% Get Netwroks names
nets = fieldnames(features);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions names
    regs = fieldnames(features.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scan names
        scans = fieldnames(features.(nets{ii}).(regs{jj}));
        % Loop over cases
        for ll = 1:size(features.(nets{ii}).(regs{jj}).(scans{1}),1)
            cVolume.(nets{ii}).(regs{jj})(ll,:) =...
                ((table2array(features.(nets{ii}).(regs{jj}).(scans{2})(ll,14)) - ...
                table2array(features.(nets{ii}).(regs{jj}).(scans{1})(ll,14)))/...
                table2array(features.(nets{ii}).(regs{jj}).(scans{1})(ll,14)))*100;
        end
    end
end
end

%% calculateSCvolume
% This function is to calculate the significance of the volume changes
%
% Input:
%   1. cVolume: The change in volumes. (structure)
%   2. test: the statistical test. (string)
%
% output:
%   1. sCvolume: The significance of changes in volumes. (structure)
%
function sCvolume = calculateSCvolume(cVolume,test)
% --Test between Scan 1 and Scan 2
sCvolume = table;
% Get Netwroks names
nets = fieldnames(cVolume);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(cVolume.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Do the statistical test
        for ll = 1:numel(nets)
            if contains(test,'signrank')
                [pvN.(regs{jj})(ll,ii),~] =...
                    signrank(cVolume.(nets{ll}).(regs{jj}),...
                    cVolume.(nets{ii}).(regs{jj}));
            elseif contains(test,'ranksum')
                [pvN.(regs{jj})(ll,ii),~] =...
                    ranksum(cVolume.(nets{ll}).(regs{jj}),...
                    cVolume.(nets{ii}).(regs{jj}));
            end
        end
    end
end

% Loop over regions
for jj = 1:numel(regs)
    % Correct
    inP = nonzeros(tril(pvN.(regs{jj}),-1));
    [~, ~, ~, sCvolume.(regs{jj})] =...
        fdr_bh(inP,0.05,'pdep','no');
end

% Add raw names
netsCN = cell(numel(netsN),numel(netsN));
for ii = 1:numel(netsN)
    for jj = ii:numel(netsN)-1
        netsCN{jj,ii} = [netsN{ii} ' Vs ' netsN{jj+1}];
    end
end
netsCN = netsCN(~cellfun('isempty',netsCN));
netsCN = netsCN(:);

sCvolume.Properties.RowNames = netsCN;
end

%% excludedIdx
% This function is to Get Idx of the excluded cases based on the quality
% control score.
%
% Input:
%   1. scores: The calculated performance scores for the cases. (structure)
%
% output:
%   1. eIdx: The excluded cases Idxs. (structure)
%
function eIdx = excludedIdx(scores)
% Get Netwroks names
nets = fieldnames(scores);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Loop over scans
    scans = fieldnames(scores.(nets{ii}).WP);
    idx = cell(1,numel(scans));
    for kk = 1:numel(scans)
        qs = scores.(nets{ii}).WP.(scans{kk}).qualityScore;
        idx{kk} = find(qs<85);
    end
    eIdx.(nets{ii}) = unique(vertcat(idx{:}));
end
end

%% featuresAE
% This function is to get features structure after excluding cases
%
% Input:
%   1. features: The extracted shape features. (structure)
%   2. eIdx: The excluded cases Idxs. (structure)
%
% output:
%   1. featuresA: features after excluding cases. (structure)
%
function featuresA = featuresAE(features,eIdx)
% Copy the Manual featres as they are
featuresA.Manual = features.Manual;
% Get Netwroks names
nets = fieldnames(eIdx);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions names
    regs = fieldnames(features.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over scans
        scans = fieldnames(features.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            featuresA.(nets{ii}).(regs{jj}).(scans{kk}) =...
                features.(nets{ii}).(regs{jj}).(scans{kk});
            featuresA.(nets{ii}).(regs{jj}).(scans{kk})(eIdx.(nets{ii}),:) = [];
        end
    end
end
end

%% iccRN
% This function is to get ICC structure after excluding cases 1000 times
% randomly
%
% Input:
%   1. features: The extracted shape features before excluding. (structure)
%   2. eIdx: The excluded cases from by quality score Idxs. (structure)
%
% output:
%   1. icc: The icc results after the random selection. (structure)
%
function icc = iccRN(features,eIdx)
% Set a seed
rng(1)
% Loop 1000 time
for rr = 1:1000
    % Get Netwroks names
    nets = fieldnames(eIdx);
    % Loop over Netwroks
    for ii = 1:numel(nets)
        % Generate random numbers to exclude
        eI.(nets{ii}) = randperm(size(features.(nets{ii}).WP.Scan1,1),...
            numel(eIdx.(nets{ii})));
    end
    % Generate new features structure after excludion
    featuresA = featuresAE(features,eI);
    % Calculate and Assign ICC after the exclusion
    icc.(['r' num2str(rr)]) = calculateICC(featuresA);
end
end

%% getScoresAE
% This function is to get the perfromance scores after excluding cases
%
% Input:
%   1. features: The extracted shape features. (structure)
%   2. scores: The calculated performance scores for the cases. (structure)
%
% output:
%   1. scoresAE: The excluded cases Idxs. (structure)
%
function scoresAE = getScoresAE(scores,eIdx)
% Get Netwroks names
nets = fieldnames(scores);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions names
    regs = fieldnames(scores.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over scans
        scans = fieldnames(scores.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            % Get into cells or structures and excluded cases
            fields = fieldnames(scores.(nets{ii}).(regs{jj}).(scans{kk}));
            for ff = 1:numel(fields)
                cf = scores.(nets{ii}).(regs{jj}).(scans{kk}).(fields{ff});
                if size(cf,1) < 2
                    % Go to the subfields first
                    subFields = fieldnames(cf);
                    for sf = 1:numel(subFields)
                        scoresAE.(nets{ii}).(regs{jj}).(scans{kk}).(fields{ff}).(subFields{sf}) =...
                            cf.(subFields{sf});
                        scoresAE.(nets{ii}).(regs{jj}).(scans{kk}).(fields{ff}).(subFields{sf})(eIdx.(nets{ii})) = [];
                    end
                else
                    scoresAE.(nets{ii}).(regs{jj}).(scans{kk}).(fields{ff}) = cf;
                    scoresAE.(nets{ii}).(regs{jj}).(scans{kk}).(fields{ff})(eIdx.(nets{ii})) = [];
                end
            end
        end
    end
end
end

%% calculateSiccBA
% This function is to calculate the significance of changes between the
% calculated ICC before and after excluding cases
%
% Input:
%   1. ICC: The calculated ICC results before excludion. (structure)
%   2. ICCAE: The calculated ICC results after excludion. (structure)
%
% output:
%   1. sICC: The significance of ICC changes. (structure)
%
function sICCBA = calculateSiccBA(ICC,ICCAE)
% --Test between Networks 95CI per feature
% Loop over Netwroks
% Get Netwroks names
nets = fieldnames(ICC);
for ii = 1:numel(nets)
    % Make table to fill later
    sICCBA.CI.(nets{ii}) = table;
    % Loop over regions
    regs = ICC.(nets{ii}).ICC.Properties.VariableNames;
    for jj = 1:numel(regs)
        % Loop over features
        fNames = ICC.(nets{ii}).ICCCI.Properties.RowNames;
        for ff = 1:numel(fNames)
            % Do the statistical test
            sICCBA.CI.(nets{ii}).(regs{jj})(ff) =...
                isempty(range_intersection(ICC.(nets{ii}).ICCCI.(regs{jj}){ff},...
                ICCAE.(nets{ii}).ICCCI.(regs{jj}){ff}));
        end
    end
    sICCBA.CI.(nets{ii}).Properties.RowNames = fNames;
end

% --Test between Networks ICC per feature
% Loop over Netwroks
% Get Netwroks names
nets = fieldnames(ICC);
for ii = 1:numel(nets)
    % Make table to fill later
    sICCBA.ICC.(nets{ii}) = table;
    % Loop over regions
    regs = ICC.(nets{ii}).ICC.Properties.VariableNames;
    for jj = 1:numel(regs)
        % Loop over features
        fNames = ICC.(nets{ii}).ICCCI.Properties.RowNames;
        for ff = 1:numel(fNames)
            % Do the statistical test
            sICCBA.ICC.(nets{ii}).(regs{jj})(ff) =...
                ICC.(nets{ii}).ICC.(regs{jj})(ff)<...
                ICCAE.(nets{ii}).ICC.(regs{jj})(ff);
        end
    end
    sICCBA.ICC.(nets{ii}).Properties.RowNames = fNames;
end
end

%% calculateSiccBAR
% This function is to calculate the significance of changes between the
% calculated ICC before (after excluding randomly 1000 times) and
% after excluding cases
%
% Input:
%   1. ICC: The icc results after the random selection. (structure)
%   2. ICCAE: The calculated ICC results after excludion. (structure)
%
% output:
%   1. sICC: The significance of ICC changes. (structure)
%
function sICC = calculateSiccBAR(ICC,ICCAE)
% --Test between Networks ICC per feature
% Loop over Netwroks
% Get Netwroks names
nets = fieldnames(ICCAE);
for ii = 1:numel(nets)
    % Make table to fill later
    sICC.(nets{ii}) = table;
    % Loop over regions
    regs = ICCAE.(nets{ii}).ICC.Properties.VariableNames;
    for jj = 1:numel(regs)
        % Loop over features
        fNames = ICCAE.(nets{ii}).ICCCI.Properties.RowNames;
        for ff = 1:numel(fNames)
            idx = zeros(1000,1);
            for rr = 1:1000
                idx(rr,:) = ICC.(['r' num2str(rr)]).(nets{ii}).ICC.(regs{jj})(ff)...
                    >=ICCAE.(nets{ii}).ICC.(regs{jj})(ff);
            end
            % Do the statistical test
            sICC.(nets{ii}).(regs{jj})(ff) = sum(idx)<50;
        end
    end
    sICC.(nets{ii}).Properties.RowNames = fNames;
end
end

%% clalcMSscores
% This function is to calculate median and Standarad deviation for multiple
% groups. For performance scores
%
% Input:
%   1. scores: the performance scores. (struct)
%
% Output:
%   1.MS: the median and standard deviation. (struct)
%
function MS = clalcMSscores(scores)
% Get Netwroks names
nets = fieldnames(scores);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(scores.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scan names
        scans = fieldnames(scores.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            MS.med.(nets{ii}).(regs{jj}).(scans{kk}) =...
                round(median(scores.(nets{ii}).(regs{jj}).(scans{kk}).WP.DSC),3);
            MS.sd.(nets{ii}).(regs{jj}).(scans{kk}) =...
                round(std(scores.(nets{ii}).(regs{jj}).(scans{kk}).WP.DSC),3);
        end
    end
end
end

%% clalcMSslices
% This function is to calculate median and Standarad deviation for multiple
% groups. For number of the included slices
%
% Input:
%   1. slices: the segmentation slices numbers. (struct)
%
% Output:
%   1.MS: the median and standard deviation. (struct)
%
function MS = clalcMSslices(slices)
% Get Netwroks names
nets = fieldnames(slices);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(slices.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scan names
        scans = fieldnames(slices.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            MS.med.(nets{ii}).(regs{jj}).(scans{kk}) =...
                round(median(slices.(nets{ii}).(regs{jj}).(scans{kk})),2);
            MS.sd.(nets{ii}).(regs{jj}).(scans{kk}) =...
                round(std(slices.(nets{ii}).(regs{jj}).(scans{kk})),2);
        end
    end
end
end

%% clalcMScvolume
% This function is to calculate median and Standarad deviation for multiple
% groups. For change in volume.
%
% Input:
%   1. cv: the chnage in volume. (struct)
%
% Output:
%   1.MS: the median and standard deviation. (struct)
%
function MS = clalcMScvolume(cv)
% Get Netwroks names
nets = fieldnames(cv);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(cv.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        
        MS.med.(nets{ii}).(regs{jj}) =...
            round(median(cv.(nets{ii}).(regs{jj})),2);
        MS.sd.(nets{ii}).(regs{jj}) =...
            round(std(cv.(nets{ii}).(regs{jj})),2);
    end
end
end

%% clalcDMSscores
% This function is to calculate difference median and Standarad deviation
%for multiple groups. For performance scores.
%
% Input:
%   1. scores: the performance scores. (struct)
%
% Output:
%   1.MS: the median and standard deviation. (struct)
%
function MS = clalcDMSscores(scores)
% Get Netwroks names
nets = fieldnames(scores);
% Loop over Netwroks
netsN = cell(numel(nets),1);
for ii = 1:numel(nets)
    % Get nets names
    netsN{ii} = nets{ii};
    % Get regions names
    regs = fieldnames(scores.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Get scan names
        scans = fieldnames(scores.(nets{ii}).(regs{jj}));
        % Calculate the difference
        dif = ((scores.(nets{ii}).(regs{jj}).(scans{end}).WP.DSC -....
            scores.(nets{ii}).(regs{jj}).(scans{1}).WP.DSC)./...
            scores.(nets{ii}).(regs{jj}).(scans{end}).WP.DSC)*100;
        % Calculate median and std
        MS.med.(nets{ii}).(regs{jj}) =...
            round(median(dif),2);
        MS.sd.(nets{ii}).(regs{jj}) =...
            round(std(dif),2);
    end
end
end
