%% ---- Prepare new manual masks ---- %%
% This function is to prepare new manual masks we are going to use in the
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
function prepNMS(masterPath)
%% Paths
paths.train.main = fullfile(masterPath,'Data','Train','Segmentation','Manual');
paths.train.original = fullfile(paths.train.main,'Original');
paths.train.adjusted = fullfile(paths.train.main,'Adjusted');
paths.test.main = fullfile(masterPath,'Data','Segmentation','Manual');
paths.test.original.Scan1 = fullfile(paths.test.main,'Original','Scan1');
paths.test.original.Scan2 = fullfile(paths.test.main,'Original','Scan2');
paths.test.adjusted.Scan1 = fullfile(paths.test.main,'Adjusted','Scan1');
paths.test.adjusted.Scan2 = fullfile(paths.test.main,'Adjusted','Scan2');
paths.clinical = fullfile(masterPath,'SheetTableAll.mat');
paths.allVois = fullfile(masterPath,'vois');
%% Load
load(paths.clinical)
%% Options
options = {'Train','Scan1','Scan2'};
%% Generate new masks
% Loop over the options
for ii = 1:numel(options)
    % Get the folder paths for original and ajusted cases
    if ii == 1
        tpo = paths.train.original;
        tpa = paths.train.adjusted;
        tpWP = fullfile(paths.train.main,'Final','WP');
        tpPZ = fullfile(paths.train.main,'Final','PZ');
        tpTZ = fullfile(paths.train.main,'Final','TZ');
        tpWPh = fullfile(paths.train.main,'Final','WPh');
        tpPZh = fullfile(paths.train.main,'Final','PZh');
        tpTZh = fullfile(paths.train.main,'Final','TZh');
        tpWPl = fullfile(paths.train.main,'Final','WPl');
        tpPZl = fullfile(paths.train.main,'Final','PZl');
        tpTZl = fullfile(paths.train.main,'Final','TZl');
        tpPZTZ = fullfile(paths.train.main,'Final','PZTZ');
        if ~exist(tpWP,'dir')
            mkdir(tpWP)
            mkdir(tpPZ)
            mkdir(tpTZ)
            mkdir(tpWPh)
            mkdir(tpPZh)
            mkdir(tpTZh)
            mkdir(tpWPl)
            mkdir(tpPZl)
            mkdir(tpTZl)
            mkdir(tpPZTZ)
        end
    else
        tpo = paths.test.original.(options{ii});
        tpa = paths.test.adjusted.(options{ii});
        tpWP = fullfile(paths.test.main,'Final',options{ii},'WP');
        tpPZ = fullfile(paths.test.main,'Final',options{ii},'PZ');
        tpTZ = fullfile(paths.test.main,'Final',options{ii},'TZ');
        tpWPh = fullfile(paths.test.main,'Final',options{ii},'WPh');
        tpPZh = fullfile(paths.test.main,'Final',options{ii},'PZh');
        tpTZh = fullfile(paths.test.main,'Final',options{ii},'TZh');
        tpWPl = fullfile(paths.test.main,'Final',options{ii},'WPl');
        tpPZl = fullfile(paths.test.main,'Final',options{ii},'PZl');
        tpTZl = fullfile(paths.test.main,'Final',options{ii},'TZl');
        tpPZTZ = fullfile(paths.test.main,'Final',options{ii},'PZTZ');
        if ~exist(tpWP,'dir')
            mkdir(tpWP)
            mkdir(tpPZ)
            mkdir(tpTZ)
            mkdir(tpWPh)
            mkdir(tpPZh)
            mkdir(tpTZh)
            mkdir(tpWPl)
            mkdir(tpPZl)
            mkdir(tpTZl)
            mkdir(tpPZTZ)
        end
    end
    % List the avilable cases in the folder
    dtpo = dir(fullfile(tpo,'*.mhd'));
    dtpa = dir(fullfile(tpa,'*.mhd'));
    % Loop over the cases inside the adjusted cases folder
    for jj = 1:numel(dtpa)
        % Display the current case name
        disp([options{ii} ' -> ' dtpa(jj).name(1:18)])
        % Current case
        FilenameO = fullfile(dtpo(jj).folder,dtpa(jj).name);
        FilenameA = fullfile(dtpa(jj).folder,dtpa(jj).name);
        currentCase = dtpa(jj).name(1:18);
        currentCaseIdx = find(contains(SheetTableAll.Sheet_ID,currentCase));
        % Read the MetaIO in matlab
        [StrDataxO, ~, ~] = elxMetaIOFileToStrDatax(FilenameO, 0);
        [StrDataxA, ~, ~] = elxMetaIOFileToStrDatax(FilenameA, 0);
        im3dO = elxStrDataxToIm3d(StrDataxO);
        im3dA = elxStrDataxToIm3d(StrDataxA);
        % Lesions index
        lesionsIdx = im3dO.Data>2;
        lesion1Idx = im3dO.Data==3;
        lesion2Idx = im3dO.Data==4;
        lesion3Idx = im3dO.Data==5;
        %% 1. Whole prostate mask including the lesions
        Prostate = im3dA;
        Prostate.Data = Prostate.Data>0;
        %% 2. Whole prostate mask excluding the lesions
        HealthyProstate = Prostate;
        HealthyProstate.Data(lesionsIdx) = 0;
        %% 3. TZ mask including the lesions
        TZ = im3dA;
        TZ.Data = TZ.Data==2;
        %% 4. TZ mask excluding the lesions
        TZHealthy = TZ;
        TZHealthy.Data(lesionsIdx) = 0;
        %% 5. PZ mask including the lesions
        PZ = im3dA;
        PZ.Data = PZ.Data==1;
        %% 6. TZ mask excluding the lesions
        PZHealthy = PZ;
        PZHealthy.Data(lesionsIdx) = 0;
        %% 7. PZ and TZ regions excluding the lesions
        PZTZ = im3dA;
        % PZTZ index
        pztzIdx = or(PZTZ.Data==1,PZTZ.Data==2);
        PZTZ.Data(~pztzIdx) = 0;
        %% 8. The lesions of the whole prostate mask
        ProstateLesions = im3dA;
        ProstateLesions.Data(:) = 0;
        ProstateLesions.Data(lesionsIdx) = 1;
        %% 9. The lesions of the TZ mask
        TZlesionsIdx = xor(TZ.Data,TZHealthy.Data);
        TZLesions = im3dA;
        TZLesions.Data(:) = 0;
        TZLesions.Data(TZlesionsIdx) = 1;
        %% 10. The lesions of the PZ mask
        PZlesionsIdx = xor(PZ.Data,PZHealthy.Data);
        PZLesions = im3dA;
        PZLesions.Data(:) = 0;
        PZLesions.Data(PZlesionsIdx) = 1;
        %% 11. The lesions indices of the whole prostate mask
        if sum(lesion1Idx,'all')>0
            ProstateLesionIdx1 = im3dA;
            ProstateLesionIdx1.Data(:) = 0;
            ProstateLesionIdx1.Data(lesion1Idx) = 1;
        elseif sum(lesion2Idx,'all')>0
            ProstateLesionIdx2 = im3dA;
            ProstateLesionIdx1.Data(:) = 0;
            ProstateLesionIdx1.Data(lesion2Idx) = 1;
        elseif sum(lesion3Idx,'all')>0
            ProstateLesionIdx3 = im3dA;
            ProstateLesionIdx3.Data(:) = 0;
            ProstateLesionIdx3.Data(lesion3Idx) = 1;
        end
        %% Biopsy data
        if exist('ProstateLesionIdx1','var')
            % If the case in the sheets list get the info
            if numel(currentCaseIdx)>0
                % Type of the biopsy
                if ~isempty(SheetTableAll.type_biopsi{currentCaseIdx})
                    BiopsyType = SheetTableAll.type_biopsi{currentCaseIdx};
                else
                    BiopsyType = 'Not avilable';
                end
                % Number of biopsies
                if ~isempty(SheetTableAll.antall_biopsier{currentCaseIdx})
                    BiopsiesNumber = SheetTableAll.antall_biopsier{currentCaseIdx};
                else
                    BiopsiesNumber ='Not avilable';
                end
                % True postive Vs. False Positive
                if ~isempty(SheetTableAll.biopsisvar{currentCaseIdx})
                    BiopsyAnswer = SheetTableAll.biopsisvar{currentCaseIdx};
                else
                    BiopsyAnswer = 'Not avilable';
                end
                % Gleason score
                if ~isempty(SheetTableAll.Dersom_positiv_biopsi_Angitt_Gleason_score_for_biopsirunden{currentCaseIdx})
                    GleasonScore = SheetTableAll.Dersom_positiv_biopsi_Angitt_Gleason_score_for_biopsirunden{currentCaseIdx};
                else
                    GleasonScore = 'Not avilable';
                end
                % Group grade
                if ~isempty(SheetTableAll.Dersom_positiv_biopsi_Angitt_Gradgruppe_for_biopsirunden{currentCaseIdx})
                    GradeGroup = SheetTableAll.Dersom_positiv_biopsi_Angitt_Gradgruppe_for_biopsirunden{currentCaseIdx};
                else
                    GradeGroup = 'Not avilable';
                end
            else
                BiopsyType = 'Not avilable';
                BiopsiesNumber ='Not avilable';
                BiopsyAnswer = 'Not avilable';
                GleasonScore = 'Not avilable';
                GradeGroup = 'Not avilable';
            end
        end
        %% Assign to a structure
        if sum(Prostate.Data,'all')>0
            vois.Prostate = Prostate;
        end
        if sum(HealthyProstate.Data,'all')>0
            vois.HealthyProstate = HealthyProstate;
        end
        if sum(TZ.Data,'all')>0
            vois.TZ = TZ;
        end
        if sum(TZHealthy.Data,'all')>0
            vois.TZHealthy = TZHealthy;
        end
        if sum(PZ.Data,'all')>0
            vois.PZ = PZ;
        end
        if sum(PZHealthy.Data,'all')>0
            vois.PZHealthy = PZHealthy;
        end
        if sum(PZTZ.Data,'all')>0
            vois.PZTZ = PZTZ;
        end
        if sum(ProstateLesions.Data,'all')>0
            vois.ProstateLesions = ProstateLesions;
        end
        if sum(TZLesions.Data,'all')>0
            vois.TZLesions = TZLesions;
        end
        if sum(PZLesions.Data,'all')>0
            vois.PZLesions = PZLesions;
        end
        if exist('ProstateLesionIdx1','var')
            vois.ProstateLesionIdx1 = ProstateLesionIdx1;
        end
        if exist('ProstateLesionIdx2','var')
            vois.ProstateLesionIdx2 = ProstateLesionIdx2;
        end
        if exist('ProstateLesionIdx3','var')
            vois.ProstateLesionIdx3 = ProstateLesionIdx3;
        end
        if exist('BiopsyType','var')
            vois.BiopsyType = BiopsyType;
        end
        if exist('BiopsiesNumber','var')
            vois.BiopsiesNumber = BiopsiesNumber;
        end
        if exist('BiopsyAnswer','var')
            vois.BiopsyAnswer = BiopsyAnswer;
        end
        if exist('GleasonScore','var')
            vois.GleasonScore = GleasonScore;
        end
        if exist('GradeGroup','var')
            vois.GradeGroup = GradeGroup;
        end
        %% Save
        % Save mat
        if ~exist(paths.allVois,'dir')
            mkdir(paths.allVois)
        end
        saveName = fullfile(paths.allVois,[currentCase '_vois.mat']);
        save(saveName,'vois')
        % Save mhd
        % names
        savemhdNameWP = fullfile(tpWP,dtpa(jj).name);
        savemhdNamePZ = fullfile(tpPZ,dtpa(jj).name);
        savemhdNameTZ = fullfile(tpTZ,dtpa(jj).name);
        savemhdNameWPh = fullfile(tpWPh,dtpa(jj).name);
        savemhdNamePZh = fullfile(tpPZh,dtpa(jj).name);
        savemhdNameTZh = fullfile(tpTZh,dtpa(jj).name);
        savemhdNamePZTZ = fullfile(tpPZTZ,dtpa(jj).name);
        savemhdNameWPl = fullfile(tpWPl,dtpa(jj).name);
        savemhdNamePZl = fullfile(tpPZl,dtpa(jj).name);
        savemhdNameTZl = fullfile(tpTZl,dtpa(jj).name);
        % from im3d to struct
        StrDataxWP = elxIm3dToStrDatax(Prostate);
        StrDataxPZ = elxIm3dToStrDatax(PZ);
        StrDataxTZ = elxIm3dToStrDatax(TZ);
        StrDataxWPh = elxIm3dToStrDatax(HealthyProstate);
        StrDataxPZh = elxIm3dToStrDatax(PZHealthy);
        StrDataxTZh = elxIm3dToStrDatax(TZHealthy);
        StrDataxPZTZ = elxIm3dToStrDatax(PZTZ);
        StrDataxWPl = elxIm3dToStrDatax(ProstateLesions);
        StrDataxPZl = elxIm3dToStrDatax(PZLesions);
        StrDataxTZl = elxIm3dToStrDatax(TZLesions);
        % from struct to MetaIO
        elxStrDataxToMetaIOFile(StrDataxWP, savemhdNameWP, 0);
        elxStrDataxToMetaIOFile(StrDataxPZ, savemhdNamePZ, 0);
        elxStrDataxToMetaIOFile(StrDataxTZ, savemhdNameTZ, 0);
        elxStrDataxToMetaIOFile(StrDataxWPh, savemhdNameWPh, 0);
        elxStrDataxToMetaIOFile(StrDataxPZh, savemhdNamePZh, 0);
        elxStrDataxToMetaIOFile(StrDataxTZh, savemhdNameTZh, 0);
        elxStrDataxToMetaIOFile(StrDataxPZTZ, savemhdNamePZTZ, 0);
        elxStrDataxToMetaIOFile(StrDataxWPl, savemhdNameWPl, 0);
        elxStrDataxToMetaIOFile(StrDataxPZl, savemhdNamePZl, 0);
        elxStrDataxToMetaIOFile(StrDataxTZl, savemhdNameTZl, 0);
        %% clean
        clearvars -except masterPath paths options SheetTableAll ii jj tpa tpo dtpa dtpo...
            tpWP tpWPh tpWPl tpPZ tpPZh tpPZl tpTZ tpTZh tpTZl tpPZTZ
    end
end