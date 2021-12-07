%% ---- Report Results ---- %%
% This function is to report the results for the
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
% 09.Nov.2020
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%   2. stats: The statistical analysis results. (structure)
%
% output:
%   Graphs, figures and tables saved to "Report" folder
%
function reportResults(masterPath,stats)
% %---- Main investigation
%% Make report folder
rpP = fullfile(masterPath,'Report','Main');
if ~exist(rpP,'dir')
    mkdir(rpP)
end
%% Scanning protocol parameters of the investigation set
protocolSummary = repProtocol(stats.protocol,rpP);
%% Scanning protocol parameters of the Training set
repProtocolTtr(stats.protocolTr,rpP);
%% Segmentation example
repSeg(masterPath,rpP)
%% ICC
repICC(stats.ICC,stats.sICC,rpP);
%% CV
repCV(stats.CV,stats.sCV,rpP);
%% Perfroamcne score
repPS(stats.pScores.pScores,stats.sScores,rpP);
%% Number of slices
repNS(stats.nSlices,stats.sSlices,rpP);
%% Volume
repVL(stats.Volume,stats.sVolume,...
    stats.cVolume,stats.sCvolume,stats.CV,protocolSummary,rpP);
%% Bland-altman analysis of Volume
blandaltmanVolume(stats.Volume,rpP);

%% Investigate Parameters
invp(stats.protocol,stats.CV);

%---- After Impleminting QC system
%% Make report folder
rpP = fullfile(masterPath,'Report','QC');
if ~exist(rpP,'dir')
    mkdir(rpP)
end
%% ICC
repICC(stats.ICCAE,stats.sICCAE,rpP);
%% CV
repCV(stats.CVAE,stats.sCVAE,rpP);
%% Perfroamcne score
repPS(stats.pScoresAE,stats.sScoresAE,rpP);
%% Number of slices
repNS(stats.nSlicesAE,stats.sSlicesAE,rpP);

%---- Before including post-processing step
%% Make report folder
rpP = fullfile(masterPath,'Report','PP');
if ~exist(rpP,'dir')
    mkdir(rpP)
end
%% ICC
repICC(stats.ICCBPP,stats.sICCBPP,rpP);
%% CV
repCV(stats.CVBPP,stats.sCVBPP,rpP);
%% Perfroamcne score
repPS(stats.pScoresBPP.pScores,stats.sScoresBPP,rpP);
%% Number of slices
repNS(stats.nSlicesBPP,stats.sSlicesBPP,rpP);
%% Volume
repVL(stats.VolumeBPP,stats.sVolumeBPP,...
    stats.cVolumeBPP,stats.sCvolumeBPP,stats.CVBPP,protocolSummary,rpP);
end

%% repProtocol
% This function is to report the scanning protocol parameters
% for Scan 1 and Scan 2 cases
%
% Input:
%   1. protocol: The scanning protocol parameters. (structure)
%   2. rpP: The report folder path. (string)
%
% Output:
%   1. protocolSummary: The summary of the scanning protocol parameters. (structure)
%
function protocolSummary = repProtocol(protocol,rpP)
% --Summary
% Make table to fill
rp.summary = table;
% Loop over scans
scans = fieldnames(protocol);
for kk = 1:numel(scans)
    % Loop over the parameters
    pars = protocol.(scans{kk}).Properties.VariableNames;
    for pp = 1:numel(pars)
        % Don't summaraize the exceptions
        exception = {'Name','ScanDate','PatientAge'};
        if ~ismember(pars{pp},exception)
            % Treat cell arrays based on uniqueness
            if iscell(protocol.(scans{kk}).(pars{pp}))
                rp.summary.(pars{pp})(kk,:) =...
                    {unique(cell2mat(protocol.(scans{kk}).(pars{pp})),'rows')};
                % Treat double arrays based on max and min
            else
                minTmp = min(protocol.(scans{kk}).(pars{pp}));
                maxTmp = max(protocol.(scans{kk}).(pars{pp}));
                % Take range if the values are different
                if minTmp ~= maxTmp
                    rp.summary.(pars{pp})(kk,:) =...
                        {[num2str(minTmp) ' - ' num2str(maxTmp)]};
                else
                    rp.summary.(pars{pp})(kk,:) = {num2str(minTmp)};
                end
            end
        end
    end
    % Assign the table rows
    rp.summary.Properties.RowNames{kk} = scans{kk};
end
%-- Patients
% Compine the scans patients
PatientAge = [protocol.Scan1.PatientAge;protocol.Scan2.PatientAge];
rp.PatientAge = ['Median = ' num2str(median(PatientAge))...
    ' Range = ' [num2str(min(PatientAge)) ' - ' num2str(max(PatientAge))]];
%-- Time interval
rp.TimeInterval = datenum(protocol.Scan2.ScanDate)-...
    datenum(protocol.Scan1.ScanDate);
rp.TimeIntervalS = ['Median = ' num2str(median(rp.TimeInterval))...
    ' Range = ' [num2str(min(rp.TimeInterval)) ' - ' num2str(max(rp.TimeInterval))]];

% Save
protocolSummary = rp;
save(fullfile(rpP,'protocolSummary.mat'),'protocolSummary')
end

%% repProtocolTr
% This function is to report the scanning protocol parameters
% for the traininng set
%
% Input:
%   1. protocol: The scanning protocol parameters. (table)
%   2. rpP: The report folder path. (string)
%
% Output:
%   protocolSummaryTr (table) saved in rpP
%
function repProtocolTtr(protocol,rpP)
% --Summary
% Make table to fill
rp.summary = table;
% Loop over the parameters
pars = protocol.Properties.VariableNames;
for pp = 1:numel(pars)
    % Don't summaraize the exceptions
    exception = {'Name','ScanDate','PatientAge'};
    if ~ismember(pars{pp},exception)
        % Treat cell arrays based on uniqueness
        if iscell(protocol.(pars{pp}))
            rp.summary.(pars{pp}) =...
                {unique(cell2mat(protocol.(pars{pp})),'rows')};
            % Treat double arrays based on max and min
        else
            minTmp = min(protocol.(pars{pp}));
            maxTmp = max(protocol.(pars{pp}));
            % Take range if the values are different
            if minTmp ~= maxTmp
                rp.summary.(pars{pp}) =...
                    {[num2str(minTmp) ' - ' num2str(maxTmp)]};
            else
                rp.summary.(pars{pp}) = {num2str(minTmp)};
            end
        end
    end
end

%-- Patients
% Compine the scans patients
PatientAge = protocol.PatientAge;
rp.PatientAge = ['Median = ' num2str(median(PatientAge))...
    ' Range = ' [num2str(min(PatientAge)) ' - ' num2str(max(PatientAge))]];

% Save
protocolSummaryTr = rp;
save(fullfile(rpP,'protocolSummaryTr.mat'),'protocolSummaryTr')
end

%% repSeg
% This function is to show examples on the manual and automated
% segmentations
%
% Input:
%   1. masterPath: The path to the master analysis file. (string)
%   2. rpP: The report folder path. (string)
%
% output:
%   Figured saved to "Report" folder
%
function repSeg(masterPath,rpP)
% Paths
imP = fullfile(masterPath,'Data','Normalized');
msPath = fullfile(masterPath,'Data','Segmentation','Manual','Final');
esPath = fullfile(masterPath,'Data','Segmentation','AutomatedPP');
% Nets
esN = dir(fullfile(esPath,'*Net*'));
esN = {esN.name}.';
mN = {'Manual'};
nets = [mN;esN];

% Go for it
% Scan1
figure('Name','Scan1 example');
for ii = 1:numel(nets)
    % Paths
    im = fullfile(imP,'Scan1','Case001_Normalized.mhd');
    if contains(nets{ii},'Manual')
        p = fullfile(msPath,'Scan1','PZ','Case001_Segmentation.mhd');
        t = fullfile(msPath,'Scan1','TZ','Case001_Segmentation.mhd');
    else
        p = fullfile(esPath,nets{ii},'PZ','Scan1','Case001_Segmentation.mhd');
        t = fullfile(esPath,nets{ii},'TZ','Scan1','Case001_Segmentation.mhd');
    end
    % Read .mhd
    [si, ~, ~] = elxMetaIOFileToStrDatax(im, 0);
    si = elxStrDataxToIm3d(si);
    [sP, ~, ~] = elxMetaIOFileToStrDatax(p, 0);
    sP = elxStrDataxToIm3d(sP);
    [sT, ~, ~] = elxMetaIOFileToStrDatax(t, 0);
    sT = elxStrDataxToIm3d(sT);
    % selected slices
    num = [8,11,15];
    % for subplot
    if ii == 1
        ff = 1:3;
    elseif ii == 2
        ff = 4:6;
    elseif ii == 3
        ff = 7:9;
    elseif ii == 4
        ff = 10:12;
    end
    % subplot
    for jj = 1:numel(num)
        % the image
        im = si.Data(:,:,num(jj));
        im = im(round(size(im,1)*.3):round(size(im,1)*.7),round(size(im,2)*.3):round(size(im,2)*.7));
        im = uint8(im);
        % the masks
        pz = sP.Data(:,:,num(jj));
        pz = pz(round(size(pz,1)*.3):round(size(pz,1)*.7),round(size(pz,2)*.3):round(size(pz,2)*.7));
        tz = sT.Data(:,:,num(jj));
        tz = tz(round(size(tz,1)*.3):round(size(tz,1)*.7),round(size(tz,2)*.3):round(size(tz,2)*.7));
        if sum(sum(pz))>0 && sum(sum(tz))>1
            pz = pz;
            tz = tz;
        elseif sum(sum(pz))>0 && sum(sum(tz))<1
            pz = pz;
            tz(1) = 1;
        elseif sum(sum(tz))>0 && sum(sum(pz))<1
            tz = tz;
            pz(1) = 1;
        end
        L = (2*pz)+tz;
        % subplot
        subplot(numel(nets),numel(num),ff(jj))
        B = labeloverlay(im,L,'Colormap','lines','Transparency',0.7);
        imshow(B,[]);
        clear im pz tz L B
    end
end
% Save
savefig(fullfile(rpP,'Scan1example.fig'))

clear im p t si sP sT num L pz tz B

% Scan2
figure('Name','Scan2 example');
for ii = 1:numel(nets)
    % Paths
    im = fullfile(imP,'Scan2','Case001_Normalized.mhd');
    if contains(nets{ii},'Manual')
        p = fullfile(msPath,'Scan2','PZ','Case001_Segmentation.mhd');
        t = fullfile(msPath,'Scan2','TZ','Case001_Segmentation.mhd');
    else
        p = fullfile(esPath,nets{ii},'PZ','Scan2','Case001_Segmentation.mhd');
        t = fullfile(esPath,nets{ii},'TZ','Scan2','Case001_Segmentation.mhd');
    end
    % Read .mhd
    [si, ~, ~] = elxMetaIOFileToStrDatax(im, 0);
    si = elxStrDataxToIm3d(si);
    [sP, ~, ~] = elxMetaIOFileToStrDatax(p, 0);
    sP = elxStrDataxToIm3d(sP);
    [sT, ~, ~] = elxMetaIOFileToStrDatax(t, 0);
    sT = elxStrDataxToIm3d(sT);
    % selected slices
    num = [9,12,16];
    % for subplot
    if ii == 1
        ff = 1:3;
    elseif ii == 2
        ff = 4:6;
    elseif ii == 3
        ff = 7:9;
    elseif ii == 4
        ff = 10:12;
    end
    % subplot
    for jj = 1:numel(num)
        % the image
        im = si.Data(:,:,num(jj));
        im = im(round(size(im,1)*.3):round(size(im,1)*.7),round(size(im,2)*.3):round(size(im,2)*.7));
        im = uint8(im);
        % the masks
        pz = sP.Data(:,:,num(jj));
        pz = pz(round(size(pz,1)*.3):round(size(pz,1)*.7),round(size(pz,2)*.3):round(size(pz,2)*.7));
        tz = sT.Data(:,:,num(jj));
        tz = tz(round(size(tz,1)*.3):round(size(tz,1)*.7),round(size(tz,2)*.3):round(size(tz,2)*.7));
        if sum(sum(pz))>0 && sum(sum(tz))>1
            pz = pz;
            tz = tz;
        elseif sum(sum(pz))>0 && sum(sum(tz))<1
            pz = pz;
            tz(1) = 1;
        elseif sum(sum(tz))>0 && sum(sum(pz))<1
            tz = tz;
            pz(1) = 1;
        end
        L = (2*pz)+tz;
        % subplot
        subplot(numel(nets),numel(num),ff(jj))
        B = labeloverlay(im,L,'Colormap','lines','Transparency',0.7);
        imshow(B,[]);
        
        clear im pz tz L B
    end
end
% Save
savefig(fullfile(rpP,'Scan2example.fig'))
end

%% repICC
% This function is to report the ICC results
%
% Input:
%   1. ICC: The calculated ICC results. (structure)
%   2. sICC: The significance of ICC changes. (structure)
%   3. rpP: The report folder path. (string)
%
function repICC(ICC,sICC,rpP)
% --Overview figures by region
% Get Netwroks names
nets = fieldnames(ICC);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions and features names
    regs = ICC.(nets{ii}).ICC.Properties.VariableNames;
    fNames = ICC.(nets{ii}).ICC.Properties.RowNames;
    % Combine data in one array
    for jj = 1:numel(regs)
        data.ICC.(regs{jj})(:,ii) = ICC.(nets{ii}).ICC.(regs{jj})';
        data.ciLow.(regs{jj})(:,ii) = cell2mat(cellfun(@(x) x(1),...
            ICC.(nets{ii}).ICCCI.(regs{jj}), 'un', 0));
        data.ciHigh.(regs{jj})(:,ii) = cell2mat(cellfun(@(x) x(2),...
            ICC.(nets{ii}).ICCCI.(regs{jj}), 'un', 0));
    end
end

% Figure
figure('Name','ICC')
% Master title
% sgtitle('ICC for regions')
% Loop over regions
for jj = 1:numel(regs)
    % Plot
    subplot(numel(regs),1,jj)
    
    % Plot groups bar
    hBar = bar(1:14,data.ICC.(regs{jj}));
    % limits
    lm = [-0.15,1];
    ylim(lm)
    % rotate x ticks
    xtickangle(15)
    % Add title and labels
    title(regs{jj})
    ylabel('ICC')
    % Turn grid
    grid on
    
    hold on
    % Calculate the number of bars in each group
    nbars = size(data.ICC.(regs{jj}), 2);
    % Get the x coordinate of the bars
    t.x = [];
    for i = 1:nbars
        t.x = [t.x ; hBar(i).XEndPoints];
    end
    % Plot the errorbars
    errorbar(t.x',data.ICC.(regs{jj}),...
        abs(data.ciLow.(regs{jj})-data.ICC.(regs{jj})),...
        abs(data.ciHigh.(regs{jj})-data.ICC.(regs{jj})),...
        'k', 'linestyle', 'none')
    % Change x ticks labels
    xticklabels(fNames')
    
    % Add significany signs
    tx = t.x';
    counter = 0;
    % between nets
    for j = 1:numel(nets)
        for i = 1:numel(fNames)
            sN = sICC.Nets.(nets{j}).(regs{jj})(i);
            if sN
                counter = counter+1;
                LL = [tx(i,1) tx(i,j)];
                % add line between them
                ld = lm(2)-lm(1);
                plot([LL(1),LL(2)],[lm(1)+ld*(j/40) lm(1)+ld*(j/40)],...
                    '.-','color',[0.333,0.333,0.333])
            end
        end
    end
    
    hold off
    % Add background
    for i = 2:2:numel(fNames)
        n = (numel(fNames)/2)/numel(fNames);
        %m = get(gca,'YTick');
        %patch([i-n,i+n,i+n,i-n],[m(1),m(1),m(end),m(end)],[0.9 0.9 0.9],'FaceAlpha',0.5)
        patch([i-n,i+n,i+n,i-n],[lm(1) lm(1) lm(2) lm(2)],[0.333 0.333 0.333],...
            'FaceAlpha',0.1,'LineStyle','none')
    end
    % Add legends
    Nnets = replace(nets,'_','');
    %     legend([Nnets','below Poor','below Moderate','below Good'],...
    %         'Location','northeastoutside')
    f = get(gca,'Children');
    if counter>0
        legend([hBar,f(numel(f)-numel(Nnets)),f(round(numel(fNames)/2)+1)],...
            [Nnets','95%CI','Significant'],'Location','northeastoutside')
    else
        legend([hBar,f(numel(f)-numel(Nnets))],...
            [Nnets','95%CI'],'Location','northeastoutside')
    end
end
% Save
savefig(fullfile(rpP,'ICC.fig'))
end

%% repCV
% This function is to report the CV results
%
% Input:
%   1. CV: The calculated CV results. (structure)
%   2. sCV: The significance of CV changes. (structure)
%   3. rpP: The report folder path. (string)
%
function repCV(CV,sCV,rpP)
% --Overview figures by region
% Get Netwroks names
nets = fieldnames(CV);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions and features names
    regs = fieldnames(CV.(nets{ii}));
    fNames = CV.(nets{ii}).(regs{1}).Properties.VariableNames;
    % Combine data in one array
    for jj = 1:numel(regs)
        for ff = 1:numel(fNames)
            data.tempF = table2array(CV.(nets{ii}).(regs{jj})(:,ff));
            for nn = 1:size(CV.(nets{ii}).(regs{jj}),1)
                data.(regs{jj})(ff,ii,nn) = data.tempF(nn);
            end
        end
    end
end
% Figure
figure('Name','CV')
% Master titles
sgtitle('CV for regions')
% Loop over regions
for jj = 1:numel(regs)
    % Plot
    subplot(numel(regs),1,jj)
    
    % Data
    dd = data.(regs{jj});
    % Limits
    lm = [min(min(min(dd))),max(max(max(dd)))];
    % Plot
    h = boxplot2(dd,1:numel(fNames),lm,fNames);
    ylim(lm)
    hold on
    % Alter linestyle and color
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:numel(nets)
        structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), h);
    end
    set([h.lwhis h.uwhis], 'linestyle', '-');
    set(h.out, 'marker', '.');
    
    % Change x ticks labels
    set(gca,'XTick',1:numel(fNames),'XTickLabel',fNames)
    xtickangle(15)
    % Titles and labels
    title(regs{jj})
    ylabel('CV')
    
    % Reference lines
    yline(0.3,'--r', 'LineWidth',1.2);
    yline(0.2,'--m', 'LineWidth',1.2);
    yline(0.1,'--b', 'LineWidth',1.2);
    
    % Add significany signs
    % between nets
    counter = 0;
    for i = 1:numel(fNames)
        sN = sCV.Nets.(regs{jj}).(fNames{i});
        for j = 1:numel(nets)-1
            if sN(j)<0.05
                counter = counter+1;
                LL = [(h.box(((numel(nets)*i)-numel(nets)+1)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+1)).XData(1))/2,...
                    (h.box(((numel(nets)*i)-numel(nets)+1+j)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+1)+j).XData(1))/2];
                % add line between them
                ld = lm(2)-lm(1);
                plot([LL(1),LL(2)],[lm(2)-ld*(j/40) lm(2)-ld*(j/40)],...
                    '.-','color',[0.333,0.333,0.333])
            end
        end
    end
    
    % Add legends
    Nnets = replace(nets,'_','');
    f = get(gca,'Children');
    legend([h.box(1:numel(Nnets),1);f(1);f(counter+1:counter+3)],...
        [Nnets','{\it p} <0.05','Below excellent','Below good','Below acceptable'],...
        'Location','northeastoutside')
    hold off
end
% Save
savefig(fullfile(rpP,'CV.fig'))
end

%% repPS
% This function is to report the performance scores results
%
% Input:
%   1. pScores: The calculated performance scores for the cases. (structure)
%   2. sScores: The significance of performance scores changes. (structure)
%   3. rpP: The report folder path. (string)
%
function repPS(pScores,sScores,rpP)
% --Overview figures by region
% Get Netwroks names
nets = fieldnames(pScores);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions and features names
    regs = fieldnames(pScores.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over scans
        scans = fieldnames(pScores.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            for nn = 1:size(pScores.(nets{ii}).(regs{jj}).(scans{kk}).totalScore,1)
                data.(regs{jj}).DSC(kk,ii,nn) =...
                    pScores.(nets{ii}).(regs{jj}).(scans{kk}).WP.DSC(nn);
                data.(regs{jj}).total(kk,ii,nn) =...
                    pScores.(nets{ii}).(regs{jj}).(scans{kk}).totalScore(nn);
            end
        end
    end
end

% Figure
figure('Name','Performance Scores')
% Master titles
sgtitle('Performance Scores')
% Loop over regions
for jj = 1:numel(regs)
    ss = 0:numel(scans);
    % ---Plot DSC
    subplot(numel(regs),2,jj+ss(jj))
    % Plot boxplots
    
    % Group names
    Nnets = replace(nets,'_','');
    % Data
    dd = data.(regs{jj}).DSC;
    % Limits
    lm = [round(min(min(min(dd))),1),1];
    % Plot
    h = boxplot2(dd,1:numel(scans),lm,Nnets');
    ylim(lm)
    
    hold on
    % Alter linestyle and color
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:numel(nets)
        structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), h);
    end
    set([h.lwhis h.uwhis], 'linestyle', '-');
    set(h.out, 'marker', '.');
    
    % Change x ticks labels
    set(gca,'XTick',1:numel(nets),'XTickLabel',scans)
    
    % Titles and labels
    title(regs{jj})
    ylabel('DSC')
    
    % Add significany signs
    % bwtween scans
    sS = sScores.Scans.WP.DSC.(regs{jj});
    
    for i = 1:numel(nets)
        if sS(i)<0.05
            LL = [(h.box(i).XData(3)+h.box(i).XData(1))/2,...
                (h.box(i+numel(nets)).XData(3)+h.box(i+numel(nets)).XData(1))/2];
            % add line between them
            ld = lm(2)-lm(1);
            plot([LL(1),LL(2)],[lm(2)-ld*(i/40) lm(2)-ld*(i/40)],...
                '.-','color',[0.333,0.333,0.333])
        end
    end
    
    % Add legend
    f = get(gca,'Children');
    legend([h.box(1:numel(Nnets),1);f(1)],[Nnets','{\it p} <0.05'],...
        'Location','northeastoutside')
    hold off
    
    % ---Plot Total
    subplot(numel(regs),numel(scans),jj+ss(jj)+1)
    
    % Plot boxplots
    % Data
    dd = data.(regs{jj}).total;
    % Limits
    lm = [floor(min(min(min(dd)))),100];
    % Plot
    h = boxplot2(dd,1:numel(scans),lm,Nnets');
    ylim(lm)
    
    hold on
    % Alter linestyle and color
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:numel(nets)
        structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), h);
    end
    set([h.lwhis h.uwhis], 'linestyle', '-');
    set(h.out, 'marker', '.');
    
    % Change x ticks labels
    set(gca,'XTick',1:numel(nets),'XTickLabel',scans)
    
    % Titles and labels
    title(regs{jj})
    ylabel('Total score')
    
    % Add significany signs
    % between scans
    sS = sScores.Scans.totalScore.(regs{jj});
    for i = 1:numel(nets)
        if sS(i)<0.05
            LL = [(h.box(i).XData(3)+h.box(i).XData(1))/2,...
                (h.box(i+numel(nets)).XData(3)+h.box(i+numel(nets)).XData(1))/2];
            % add line between them
            ld = lm(2)-lm(1);
            plot([LL(1),LL(2)],[lm(2)-ld*(i/40) lm(2)-ld*(i/40)],...
                '.-','color',[0.333,0.333,0.333])
        end
    end
    % Add legend
    f = get(gca,'Children');
    legend([h.box(1:numel(Nnets),1);f(1)],[Nnets','{\it p} <0.05'],...
        'Location','northeastoutside')
    
    hold off
    
    %     % Change lines width for all graph
    %     lines = findobj(gcf,'Type','Line');
    %     for i = 1:numel(lines)
    %         lines(i).LineWidth = 2.0;
    %     end
end
% Save
savefig(fullfile(rpP,'Performance Scores.fig'))

% Figure
figure('Name','Performance DSC Scores')
% Master titles
% sgtitle('Performance Scores')
% Loop over regions
for jj = 1:numel(regs)
    % ---Plot DSC
    subplot(numel(regs),1,jj)
    % Plot boxplots
    
    % Group names
    Nnets = replace(nets,'_','');
    % Data
    dd = data.(regs{jj}).DSC;
    % Limits
    lm = [round(min(min(min(dd))),1),1.15];
    % Plot
    h = boxplot2(dd,1:numel(scans),lm,Nnets');
    ylim(lm)
    
    hold on
    % Alter linestyle and color
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:numel(nets)
        structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), h);
    end
    set([h.lwhis h.uwhis], 'linestyle', '-');
    set(h.out, 'marker', '.');
    
    % Change x ticks labels
    set(gca,'XTick',1:numel(nets),'XTickLabel',scans)
    
    % Titles and labels
    title(regs{jj})
    ylabel('DSC')
    
    % Add significany signs
    % bwtween scans
    sS = sScores.Scans.WP.DSC.(regs{jj});
    
    for i = 1:numel(nets)
        if sS(i)<0.05
            LL = [(h.box(i).XData(3)+h.box(i).XData(1))/2,...
                (h.box(i+numel(nets)).XData(3)+h.box(i+numel(nets)).XData(1))/2];
            % add line between them
            ld = lm(2)-lm(1);
            plot([LL(1),LL(2)],[lm(2)-ld*(i/40) lm(2)-ld*(i/40)],...
                '.-','color',[0.333,0.333,0.333])
        end
    end
    
    % between nets
    counter = 0;
    for i = 1:numel(scans)
        sN = sScores.(['NetsScan' num2str(i)]).WP.DSC.(regs{jj});
        for j = 1:numel(nets)-1
            if sN(j)<0.05
                counter = counter+1;
                LL = [(h.box(((numel(nets)*i)-numel(nets)+1)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+1)).XData(1))/2,...
                    (h.box(((numel(nets)*i)-numel(nets)+1+j)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+1)+j).XData(1))/2];
                % add line between them
                plot([LL(1),LL(2)],[lm(1)+ld*(j/40) lm(1)+ld*(j/40)],...
                    '.--','color',[0.333,0.333,0.333])
            end
        end
        for j = numel(nets)
            if sN(j)<0.05
                counter = counter+1;
                LL = [(h.box(((numel(nets)*i)-numel(nets)+2)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+2)).XData(1))/2,...
                    (h.box(((numel(nets)*i)-numel(nets)+j)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets))+j).XData(1))/2];
                % add line between them
                plot([LL(1),LL(2)],[lm(1)+ld*(j/40) lm(1)+ld*(j/40)],...
                    '.--','color',[0.333,0.333,0.333])
                
            end
        end
    end
    
    % Add legend
    f = get(gca,'Children');
    legend([h.box(1:numel(Nnets),1);f(1);f(counter+1)],...
        [Nnets','{\it p} <0.05 between netwroks','{\it p} <0.05 between scans'],...
        'Location','northeastoutside')
    hold off
    
    %     % Change lines width for all graph
    %     lines = findobj(gcf,'Type','Line');
    %     for i = 1:numel(lines)
    %         lines(i).LineWidth = 2.0;
    %     end
end
% Save
savefig(fullfile(rpP,'Performance DSC Scores.fig'))
end

%% repNS
% This function is to report the number of slices results
%
% Input:
%   1. nSlices: The numbers of the included slices. (structure)
%   2. sSlices: The significance of slices numbers change. (structure)
%   3. rpP: The report folder path. (string)
%
function repNS(nSlices,sSlices,rpP)
% --Overview figures by region
% Get Netwroks names
nets = fieldnames(nSlices);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions and features names
    regs = fieldnames(nSlices.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        % Loop over scans
        scans = fieldnames(nSlices.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            for nn = 1:size(nSlices.(nets{ii}).(regs{jj}).(scans{kk}),1)
                data.(regs{jj})(kk,ii,nn) =...
                    nSlices.(nets{ii}).(regs{jj}).(scans{kk})(nn);
            end
        end
    end
end

% Figure
figure('Name','Number of slices')
% Master title
sgtitle('Number of mask slices')
% Loop over regions
for jj = 1:numel(regs)
    subplot(numel(regs),1,jj)
    % Plot boxplots
    
    % Group names
    Nnets = replace(nets,'_','');
    % Data
    dd = data.(regs{jj});
    % Limits
    lm = [min(min(min(dd)))-2,max(max(max(dd)))+2];
    % Plot
    h = boxplot2(dd,1:numel(scans),lm,Nnets');
    ylim(lm)
    
    hold on
    % Alter linestyle and color
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:numel(nets)
        structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), h);
    end
    set([h.lwhis h.uwhis], 'linestyle', '-');
    set(h.out, 'marker', '.');
    
    % Change x ticks labels
    set(gca,'XTick',1:numel(nets),'XTickLabel',scans)
    
    % Titles and labels
    title(regs{jj})
    ylabel('Number of mask slices')
    
    % Add significany signs
    % bwtween scans
    sS = sSlices.Scans.(regs{jj});
    for i = 1:numel(nets)
        if sS(i)<0.05
            LL = [(h.box(i).XData(3)+h.box(i).XData(1))/2,...
                (h.box(i+numel(nets)).XData(3)+h.box(i+numel(nets)).XData(1))/2];
            % add line between them
            ld = lm(2)-lm(1);
            plot([LL(1),LL(2)],[lm(2)-ld*(i/40) lm(2)-ld*(i/40)],...
                '.-','color',[0.333,0.333,0.333])
        end
    end
    % between nets
    counter = 0;
    for i = 1:numel(scans)
        sN = sSlices.(['NetsScan' num2str(i)]).(regs{jj});
        for j = 1:numel(nets)-1
            if sN(j)<0.05
                counter = counter+1;
                LL = [(h.box(((numel(nets)*i)-numel(nets)+1)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+1)).XData(1))/2,...
                    (h.box(((numel(nets)*i)-numel(nets)+1+j)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+1)+j).XData(1))/2];
                % add line between them
                plot([LL(1),LL(2)],[lm(1)+ld*(j/40) lm(1)+ld*(j/40)],...
                    '.--','color',[0.333,0.333,0.333])
            end
        end
    end
    
    % Add legend
    f = get(gca,'Children');
    legend([h.box(1:numel(Nnets),1);f(counter+1);f(1)],...
        [Nnets','{\it p} <0.05 between scans','{\it p} <0.05 between methods'],...
        'Location','northeastoutside')
    hold off
    
    %         % Change lines width for all graph
    %         lines = findobj(gcf,'Type','Line');
    %         for i = 1:numel(lines)
    %             lines(i).LineWidth = 2.0;
    %         end
end
% Save
savefig(fullfile(rpP,'Number of slices.fig'))
end

%% repVL
% This function is to report the Volume related results
%
% Input:
%   1. Volume: The volumes. (structure)
%   2. sVolume: The significance of volumes change. (structure)
%   3. cVolume: The change in volumes. (structure)
%   4. sCvolume: The significance of changes in volumes. (structure)
%   5. CV: The calculated CV results. (structure)
%   6. protocol: The summary of the scanning protocol parameters. (structure)
%   7. rpP: The report folder path. (string)
%
function repVL(Volume,sVolume,cVolume,sCvolume,CV,protocol,rpP)
% --Overview figures by region
% Get Netwroks names
nets = fieldnames(Volume);
% Loop over Netwroks
for ii = 1:numel(nets)
    % Get regions and features names
    regs = fieldnames(Volume.(nets{ii}));
    % Loop over regions
    for jj = 1:numel(regs)
        %-- Volumes
        % Loop over scans
        scans = fieldnames(Volume.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            for nn = 1:size(Volume.(nets{ii}).(regs{jj}).(scans{kk}),1)
                data.volume.(regs{jj})(kk,ii,nn) =...
                    Volume.(nets{ii}).(regs{jj}).(scans{kk})(nn);
            end
        end
        %-- Change in volume
        data.cVolume.(regs{jj})(:,ii) =...
            cVolume.(nets{ii}).(regs{jj});
        
    end
end

%--Volumes
% Figure
figure('Name','Volumes in scans')
% Master title
sgtitle('Volumes in scans')
% Loop over regions
for jj = 1:numel(regs)
    subplot(numel(regs),1,jj)
    % Plot boxplots
    % Group names
    Nnets = replace(nets,'_','');
    % Data
    dd = data.volume.(regs{jj});
    % Limits
    lm = [min(min(min(dd)))-(1.5*10^4),max(max(max(dd)))+(1.5*10^4)];
    % Plot
    h = boxplot2(dd,1:numel(scans),lm,Nnets');
    ylim(lm)
    
    hold on
    % Alter linestyle and color
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:numel(nets)
        structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), h);
    end
    set([h.lwhis h.uwhis], 'linestyle', '-');
    set(h.out, 'marker', '.');
    
    % Change x ticks labels
    set(gca,'XTick',1:numel(nets),'XTickLabel',scans)
    
    % Titles and labels
    title(regs{jj})
    ylabel('Volume')
    
    % Add significany signs
    % bwtween scans
    sS = sVolume.Scans.(regs{jj});
    for i = 1:numel(nets)
        if sS(i)<0.05
            LL = [(h.box(i).XData(3)+h.box(i).XData(1))/2,...
                (h.box(i+numel(nets)).XData(3)+h.box(i+numel(nets)).XData(1))/2];
            % add line between them
            ld = lm(2)-lm(1);
            plot([LL(1),LL(2)],[lm(2)-ld*(i/40) lm(2)-ld*(i/40)],...
                '.-','color',[0.333,0.333,0.333])
        end
    end
    % between nets
    counter = 0;
    for i = 1:numel(scans)
        sN = sVolume.(['NetsScan' num2str(i)]).(regs{jj});
        for j = 1:numel(nets)-1
            if sN(j)<0.05
                counter = counter+1;
                LL = [(h.box(((numel(nets)*i)-numel(nets)+1)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+1)).XData(1))/2,...
                    (h.box(((numel(nets)*i)-numel(nets)+1+j)).XData(3)+...
                    h.box(((numel(nets)*i)-numel(nets)+1)+j).XData(1))/2];
                % add line between them
                plot([LL(1),LL(2)],[lm(1)+ld*(j/40) lm(1)+ld*(j/40)],...
                    '.--','color',[0.333,0.333,0.333])
            end
        end
    end
    
    % Add legend
    f = get(gca,'Children');
    legend([h.box(1:numel(Nnets),1);f(counter+1);f(1)],...
        [Nnets','{\it p} <0.05 between scans','{\it p} <0.05 between methods'],...
        'Location','northeastoutside')
    hold off
    
    %         % Change lines width for all graph
    %         lines = findobj(gcf,'Type','Line');
    %         for i = 1:numel(lines)
    %             lines(i).LineWidth = 2.0;
    %         end
end
% Save
savefig(fullfile(rpP,'Volumes in scans.fig'))

%--Change in volume between scans
% Figure
figure('Name','Change in volume between scans')
% Master title
sgtitle('Change in volume between scans')
% Loop over regions
for jj = 1:numel(regs)
    subplot(numel(regs),1,jj)
    % Plot boxplots
    % Group names
    Nnets = replace(nets,'_','');
    % Data
    dd = data.cVolume.(regs{jj});
    % Limits
    lm = [-100,max(max(max(dd)))];
    %     lm = [min(min(min(dd))),max(max(max(dd)))];
    % Plot
    %     cmap = colormap(lines(numel(nets)));
    cmap = colormap(lines(numel(nets)));
    boxplot(dd,'Labels',Nnets','ColorGroup',cmap);
    
    ylim(lm)
    grid on
    
    hold on
    % Change x ticks labels
    set(gca,'XTick',1:numel(nets),'XTickLabel',Nnets')
    
    % Titles and labels
    title(regs{jj})
    ylabel('Volume 2 - 1 (%)')
    
    % Add mean
    gxt = get(gca,'XTick');
    tx = zeros(numel(nets),1);
    for i = 1:numel(nets)
        tx(i) = gxt(i);
    end
    m = mean(dd);
    plot(tx ,m,'*k')
    
    % Add eros line
    yline(0,'--m', 'LineWidth',2);
    
    % Add significany signs
    % between nets
    counter = 0;
    sN = sCvolume.(regs{jj});
    for i = 2:numel(nets)
        if sN(i)<0.05
            counter = counter+1;
            LL = [gxt(1),gxt(i)];
            ld = lm(2)-lm(1);
            % add line between them
            plot([LL(1),LL(2)],[lm(1)+ld*(i/40) lm(1)+ld*(i/40)],...
                '.--','color',[0.333,0.333,0.333])
        end
    end
    
    % Add legend
    f = get(gca,'Children');
    legend([f(1),f(counter+1)],{'{\it p} <0.05','Zero'},...
        'Location','northeastoutside')
    
    hold off
    
    %     % Change lines width for all graph
    %     linesc = findobj(gca,'Type','Line');
    %     for i = 1:numel(linesc)
    %         linesc(i).LineWidth = 2.0;
    %     end
end
% Save
savefig(fullfile(rpP,'Change in volume between scans.fig'))

% %-- Correlation between voxel volume CV and interval scanning days
% figure('Name','Volume Vs. Days');
% % Master title
% sgtitle('Volume Vs. Days')
% % Loop over regions
% for jj = 1:numel(regs)
%     subplot(numel(regs),1,jj)
%     % Get days and sort them
%     [days,dIdx] = sort(protocol.TimeInterval);
%     % Get volumes CV
%     vCV = CV.Manual.(regs{jj}).VoxelVolume(dIdx);
%     % Get fitted values
%     [coeffs,S] = polyfit(days,vCV,1);
%     fittedX = linspace(0,max(days));
%     [fittedY,delta] = polyval(coeffs,fittedX,S);
%     % Plot the fitted line
%     scatter(days,vCV,'filled','k','MarkerFaceAlpha',0.4)
%     grid on
%     hold on;
%     plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
%     p = plot(fittedX, fittedY+2*delta,'m--',fittedX, fittedY-2*delta,'m--');
%     xlabel('Time interval between scans (Days)')
%     ylabel('Voxel volume CV')
%     title(regs{jj})
%     %     refline([1,0])
%     % Add Legend
%     set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     legend('Data','Linear Fit','95% Prediction Interval','Unity Line','Location','northeastoutside')
%     % Add rho
%     [r,p]= corr(days,vCV);
%     gxt = get(gca,'XTick');
%     gyt = get(gca,'YTick');
%     text(gxt(end-1),gyt(end),['rho = ',num2str(round(r,3)),'   p-val = ',num2str(round(p,3))])
%     hold off;
% end
% savefig(fullfile(rpP,'Volume Vs Days.fig'))
end

%% invp
% This function is to investigate the relationship between volume CV and
% different protocol parameters
%
% Input:
%   1. protocol: The summary of the scanning protocol parameters. (structure)
%   2. CV: The calculated CV results. (structure)
%
function invp(protocol,CV)
% Loop over Netwroks
nets = fieldnames(CV);
for ii = 1:numel(nets)
    % Get regions and features names
    regs = fieldnames(CV.(nets{ii}));
end
%-- Correlation between voxel volume CV and patients age
figure('Name','Volume Vs. Age');
% Master title
sgtitle('Volume Vs. Age')
% Loop over regions
for jj = 1:numel(regs)
    subplot(numel(regs),1,jj)
    % Get days and sort them
    [days,dIdx] = sort(protocol.Scan1.PatientAge);
    % Get volumes CV
    vCV = CV.Manual.(regs{jj}).VoxelVolume(dIdx);
    % Get fitted values
    [coeffs,S] = polyfit(days,vCV,1);
    fittedX = linspace(0,max(days));
    [fittedY,delta] = polyval(coeffs,fittedX,S);
    % Plot the fitted line
    scatter(days,vCV,'filled','k','MarkerFaceAlpha',0.4)
    grid on
    hold on;
    plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
    p = plot(fittedX, fittedY+2*delta,'m--',fittedX, fittedY-2*delta,'m--');
    xlabel('Patient age (years)')
    ylabel('Voxel volume CV')
    title(regs{jj})
    %     refline([1,0])
    % Add Legend
    set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('Data','Linear Fit','95% Prediction Interval','Unity Line','Location','northeastoutside')
    % Add rho
    [r,p]= corr(days,vCV);
    gxt = get(gca,'XTick');
    gyt = get(gca,'YTick');
    text(gxt(end-1),gyt(end),['rho = ',num2str(round(r,3)),'   p-val = ',num2str(round(p,3))])
    hold off;
end

%-- Correlation between voxel volume CV and number of averages
figure('Name','Volume Vs. Averages number');
% Master title
sgtitle('Volume Vs. Averages number')
% Loop over regions
for jj = 1:numel(regs)
    subplot(numel(regs),1,jj)
    % Get days and sort them
    [days,dIdx] = sort(protocol.Scan2.NumberOfAverages);
    % Get volumes CV
    vCV = CV.Manual.(regs{jj}).VoxelVolume(dIdx);
    % Get fitted values
    [coeffs,S] = polyfit(days,vCV,1);
    fittedX = linspace(0,max(days));
    [fittedY,delta] = polyval(coeffs,fittedX,S);
    % Plot the fitted line
    scatter(days,vCV,'filled','k','MarkerFaceAlpha',0.4)
    grid on
    hold on;
    plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
    p = plot(fittedX, fittedY+2*delta,'m--',fittedX, fittedY-2*delta,'m--');
    xlabel('Number of averages in Scan 2')
    ylabel('Voxel volume CV')
    title(regs{jj})
    %     refline([1,0])
    % Add Legend
    set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('Data','Linear Fit','95% Prediction Interval','Unity Line','Location','northeastoutside')
    % Add rho
    [r,p]= corr(days,vCV);
    gxt = get(gca,'XTick');
    gyt = get(gca,'YTick');
    text(gxt(end-1),gyt(end),['rho = ',num2str(round(r,3)),'   p-val = ',num2str(round(p,3))])
    hold off;
end

%-- Correlation between voxel volume CV and Repertition Time difference
figure('Name','Volume Vs. RT difference');
% Master title
sgtitle('Volume Vs. RT difference')
% Loop over regions
for jj = 1:numel(regs)
    subplot(numel(regs),1,jj)
    % Get days and sort them
    [days,dIdx] = sort(abs(protocol.Scan1.RepetitionTime - protocol.Scan2.RepetitionTime));
    % Get volumes CV
    vCV = CV.Manual.(regs{jj}).VoxelVolume(dIdx);
    % Get fitted values
    [coeffs,S] = polyfit(days,vCV,1);
    fittedX = linspace(0,max(days));
    [fittedY,delta] = polyval(coeffs,fittedX,S);
    % Plot the fitted line
    scatter(days,vCV,'filled','k','MarkerFaceAlpha',0.4)
    grid on
    hold on;
    plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
    p = plot(fittedX, fittedY+2*delta,'m--',fittedX, fittedY-2*delta,'m--');
    xlabel('absolute RT difference')
    ylabel('Voxel volume CV')
    title(regs{jj})
    %     refline([1,0])
    % Add Legend
    set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('Data','Linear Fit','95% Prediction Interval','Unity Line','Location','northeastoutside')
    % Add rho
    [r,p]= corr(days,vCV);
    gxt = get(gca,'XTick');
    gyt = get(gca,'YTick');
    text(gxt(end-1),gyt(end),['rho = ',num2str(round(r,3)),'   p-val = ',num2str(round(p,3))])
    hold off;
end

%-- Correlation between voxel volume CV and Repertition Slices number difference
figure('Name','Volume Vs. Slices number difference');
% Master title
sgtitle('Volume Vs. Slices number difference')
% Loop over regions
for jj = 1:numel(regs)
    subplot(numel(regs),1,jj)
    % Get days and sort them
    [days,dIdx] = sort(abs(protocol.Scan1.SlicesNumber - protocol.Scan2.SlicesNumber));
    % Get volumes CV
    vCV = CV.Manual.(regs{jj}).VoxelVolume(dIdx);
    % Get fitted values
    [coeffs,S] = polyfit(days,vCV,1);
    fittedX = linspace(0,max(days));
    [fittedY,delta] = polyval(coeffs,fittedX,S);
    % Plot the fitted line
    scatter(days,vCV,'filled','k','MarkerFaceAlpha',0.4)
    grid on
    hold on;
    plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
    p = plot(fittedX, fittedY+2*delta,'m--',fittedX, fittedY-2*delta,'m--');
    xlabel('Absolute slices number difference')
    ylabel('Voxel volume CV')
    title(regs{jj})
    %     refline([1,0])
    % Add Legend
    set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('Data','Linear Fit','95% Prediction Interval','Unity Line','Location','northeastoutside')
    % Add rho
    [r,p]= corr(days,vCV);
    gxt = get(gca,'XTick');
    gyt = get(gca,'YTick');
    text(gxt(end-1),gyt(end),['rho = ',num2str(round(r,3)),'   p-val = ',num2str(round(p,3))])
    hold off;
end

%-- Correlation between RT and Slices number
figure('Name','RT Vs. Slices number');
% Loop over regions
% Get days and sort them
[slices,dIdx] = sort(protocol.Scan2.SlicesNumber);
% Get RT
RT = protocol.Scan2.RepetitionTime(dIdx);
% Get fitted values
[coeffs,S] = polyfit(slices,RT,1);
fittedX = linspace(0,max(slices));
[fittedY,delta] = polyval(coeffs,fittedX,S);
% Plot the fitted line
scatter(slices,RT,'filled','k','MarkerFaceAlpha',0.4)
grid on
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
p = plot(fittedX, fittedY+2*delta,'m--',fittedX, fittedY-2*delta,'m--');
xlabel('Slices number in Scan 2')
ylabel('RT')
title('RT Vs. Slices number')
%     refline([1,0])
% Add Legend
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('Data','Linear Fit','95% Prediction Interval','Unity Line','Location','northeastoutside')
% Add rho
[r,p]= corr(slices,RT);
gxt = get(gca,'XTick');
gyt = get(gca,'YTick');
text(gxt(end-1),gyt(end),['rho = ',num2str(round(r,3)),'   p-val = ',num2str(round(p,3))])
hold off;

%-- Correlation between Number of averages and Slices number
figure('Name','Number of averages Vs. Slices number');
% Loop over regions
% Get days and sort them
[slices,dIdx] = sort(protocol.Scan2.SlicesNumber);
% Get Number of averages
NA = protocol.Scan2.NumberOfAverages(dIdx);
% Get fitted values
[coeffs,S] = polyfit(slices,NA,1);
fittedX = linspace(0,max(slices));
[fittedY,delta] = polyval(coeffs,fittedX,S);
% Plot the fitted line
scatter(slices,NA,'filled','k','MarkerFaceAlpha',0.4)
grid on
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
p = plot(fittedX, fittedY+2*delta,'m--',fittedX, fittedY-2*delta,'m--');
xlabel('Slices number in Scan 2')
ylabel('Number of averages')
title('Number of averages Vs. Slices number')
%     refline([1,0])
% Add Legend
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('Data','Linear Fit','95% Prediction Interval','Unity Line','Location','northeastoutside')
% Add rho
[r,p]= corr(slices,NA);
gxt = get(gca,'XTick');
gyt = get(gca,'YTick');
text(gxt(end-1),gyt(end),['rho = ',num2str(round(r,3)),'   p-val = ',num2str(round(p,3))])
hold off;

%-- Correlation between Number of averages and RT
figure('Name','Number of averages Vs. RT');
% Loop over regions
% Get days and sort them
[NA,dIdx] = sort(protocol.Scan2.NumberOfAverages);
% Get Number of averages
RT = protocol.Scan2.RepetitionTime(dIdx);
% Get fitted values
[coeffs,S] = polyfit(NA,RT,1);
fittedX = linspace(0,max(NA));
[fittedY,delta] = polyval(coeffs,fittedX,S);
% Plot the fitted line
scatter(NA,RT,'filled','k','MarkerFaceAlpha',0.4)
grid on
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
p = plot(fittedX, fittedY+2*delta,'m--',fittedX, fittedY-2*delta,'m--');
xlabel('Number of averages')
ylabel('RT')
title('Number of averages Vs. Slices number')
%     refline([1,0])
% Add Legend
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('Data','Linear Fit','95% Prediction Interval','Unity Line','Location','northeastoutside')
% Add rho
[r,p]= corr(NA,RT);
gxt = get(gca,'XTick');
gyt = get(gca,'YTick');
text(gxt(end-1),gyt(end),['rho = ',num2str(round(r,3)),'   p-val = ',num2str(round(p,3))])
hold off;
end

%% blandaltmanVolume
% This function is to report the Bland-Altman analysis of Volume
%
% Input:
%   1. Volume: The volumes. (structure)
%   2. rpP: The report folder path. (string)
%
function blandaltmanVolume(Volume,rpP)
% -- Scans
% Loop over Netwroks
nets = fieldnames(Volume);
for ii = 1:numel(nets)
    % Get regions
    regs = fieldnames(Volume.(nets{ii}));
    for jj = 1:numel(regs)
        % Figure
        figure('Name',['Bland-Altman of Volume ' (nets{ii}) ' ' (regs{jj})]);
        
        % Limits of agreement
        S1 = Volume.(nets{ii}).(regs{jj}).Scan1/1000; %/100 to get mL
        S2 = Volume.(nets{ii}).(regs{jj}).Scan2/1000;
        
        s = ba(S2,S1);
        muD = s.muD; % mean difference (bias)
        sD = s.sD; % standard deviation of difference
        loa = s.loa; % limits of agreement
        % Plot
        ba(gcf,S2,S1, 'XName','Scan 2', 'YName','Scan 1', ...
            'PlotMeanDifference',true, 'PlotStatistics','extended');
        title([nets{ii} ' ' (regs{jj})])
        % Save
        savefig(fullfile(rpP,'BA','scans',['Bland-Altman of Volume ' (nets{ii}) ' ' (regs{jj}) '.fig']))
    end
end
%-- Methods
% Loop over Netwroks
nets = fieldnames(Volume);
% Exclude Manual
nets = nets(~ismember(nets,'Manual'));
for ii = 1:numel(nets)
    % Get regions
    regs = fieldnames(Volume.(nets{ii}));
    for jj = 1:numel(regs)
        % Get Scans
        scans = fieldnames(Volume.(nets{ii}).(regs{jj}));
        for kk = 1:numel(scans)
            
            % Figure
            figure('Name',...
                ['Bland-Altman of Volume ' (nets{ii}) ' ' (regs{jj}) ' ' scans{kk}]);
            
            % Limits of agreement
            S1 = Volume.Manual.(regs{jj}).(scans{kk})/1000;
            S2 = Volume.(nets{ii}).(regs{jj}).(scans{kk})/1000;
            
            s = ba(S2,S1);
            muD = s.muD; % mean difference (bias)
            sD = s.sD; % standard deviation of difference
            loa = s.loa; % limits of agreement
            % Plot
            baM(gcf,S2,S1, 'XName',nets{ii}, 'YName','Manual', ...
                'PlotMeanDifference',true, 'PlotStatistics','extended');
            title([nets{ii} ' ' regs{jj} '' scans{kk}])
            % Save
            savefig(fullfile(rpP,'BA','methods',scans{kk},...
                ['Bland-Altman of Volume ' nets{ii} ' ' regs{jj} ' ' scans{kk} '.fig']))
        end
    end
end
end

