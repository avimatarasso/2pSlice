%% initialize wash
% This script is to analyze 2p slice sessions where a wash was done
%assumes there is a folder where delays are stored
clear all

addpath('F:\code\2pSliceAnalysis\oir2stdData-master');
addpath('C:\Users\avima\OneDrive\Documents\GitHub\2pSlice');
pathToFile = 'F:\2P-slice exps\NE DA proj\dlight\washes\saved ROIs';%F:\2P-slice exps\NE DA proj\GRABNE\washes\';%'F:\2P-slice exps\022222\';%'F:\2P-slice exps\022222\105304-3_S1_dL_CA1_20hz_30s_15sdelay.oir';
cd(pathToFile)
stimON    = 0; %this script is Wash only

%only need to change if you want to put manual title names.
manualTitles = 0;
TitleNames = {'30uM'};%{'100uM 15s_2','100uM 15s 0001','100uM 15s 0002'} ;%{'Start','End'};

close all 
% Can loop through fileOI 

%dose = {'100uM_15s','100uM_15s_0001','100uM_15s_0002'} 

% TO ADD: - a way to not do all files that come up in workFiles
%         - a way to add in specific delays for each session(import data from
%         csv? (handDelay)
%       -

ligand = 'DA'; %beta
doses = {'3nM', '10nM', '30nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'};
for dd =1:length(doses)
    filePattern = doses{dd};
    workFiles = dir(['*' filePattern '*.mat']);
    workFiles = workFiles(~endsWith({workFiles.name}, '_Analyzed.mat'));
    blVal     = zeros(length(workFiles),1);  
    peakVal   = zeros(length(workFiles),1);
    semarr    = peakVal;
    concNames = {}; 
    tic
    
for fileN = 1:length(workFiles)
     clearvars -except dd fileN semarr concNames peakVal blVal workFiles ...
         filePattern doses ligand pathToFile manualTitles averageFluor
     
    fileName = workFiles(fileN).name;
    splitStr = regexp(fileName,'_','split');
    len  = length(splitStr);
    doseOI = splitStr(contains(splitStr,ligand)); 
    doseOI = doseOI{1}; dose = doseOI(length(ligand)+1:end);
    slice = splitStr{2};
    date = splitStr(contains(splitStr,'22')); date = date{1};
    handDelay = 12;
    % load delay and file
    
    load([pathToFile(1) ':\2P-slice exps\NE DA proj\delays\delay' date '.mat']); %get delay
    load(fileName)
    if fileN == 1 % initialize
        averageFluor = zeros(276,length(workFiles)); %CUSTOMIZE %size(sessionImgs,3),length(workFiles));
    end
    
    laterDate = compareDates(fileName, '220214'); % '220214' is when i switched to resonant
    if laterDate
        fs = 3; %frame avging of 10hz
        artifactThreshold = 300; 
    else
        fs = 1/1.088;
        artifactThreshold = 500; 
    end
    dsFactor = 3/1.088;
    totDelay = delay + handDelay;
    totDelaySamp = ceil(totDelay*fs); %the total delay in samples

    %samples before for baseline
    sampBefore = ceil(30*fs); %30 s
    baselineSamples = totDelaySamp - sampBefore:totDelaySamp; 
    if exist('sessionImgs','var')
        meanBL = mean(sessionImgs(:,:,baselineSamples),3); meanBL = imgaussfilt(meanBL,2);
    else
        warning([fileName ' has no sessionImgs saved!!!'])
        continue
    end
 
    %find sum of fluorescence per area
    pixels = sum(ROImask(:));
    contourImgs = (sessionImgs.*ROImask)/pixels; 
    % find fluorescence change per area
    tmpAvgFluor = squeeze(sum(contourImgs,1)); tmpAvgFluor = squeeze(sum(tmpAvgFluor,1));
    
    if length(tmpAvgFluor)> 300
        sample2End = ceil(size(sessionImgs,3)/3);
        if sample2End >276
            sample2End = 276;
        tmpAvgFluor = decimate(tmpAvgFluor,3); tmpAvgFluor= tmpAvgFluor(3:sample2End+2);
        elseif sample2End < 276 && size(sessionImgs,3)>276
             tmpAvgFluor = decimate(tmpAvgFluor,3); 
             bltmp(1:276 - sample2End) =  mean(tmpAvgFluor(1:276 - sample2End),1);
             if size(tmpAvgFluor,1)==1
                 tmpAvgFluor = tmpAvgFluor';
             end
            tmpAvgFluor= [bltmp'; tmpAvgFluor(1:sample2End)];
        end
    elseif size(sessionImgs,3)<276
         if size(tmpAvgFluor,1)==1
                 tmpAvgFluor = tmpAvgFluor';
         end
        bltmp(1:276-length(tmpAvgFluor)) =  mean(tmpAvgFluor(1:276-length(tmpAvgFluor)),1);
        tmpAvgFluor= [bltmp'; tmpAvgFluor];            
    end
    
    tmpAvgFluorNoSpikes = medfilt1(tmpAvgFluor',3);
    f0  = mean(tmpAvgFluorNoSpikes(baselineSamples),1); f0 = mean((f0(:))~=0);
    
    averageFluor(:,fileN) = tmpAvgFluorNoSpikes;
    
    % Calculate df/f
    dff = (tmpAvgFluorNoSpikes - f0)./f0;
    if size(dff,1) == 1
        dff = dff';
    end
    %[max_dff, max] = max(dff);
    tmpBL = mean(dff(baselineSamples,:));
    blVal(fileN) =tmpBL;
    % Create a time array
    %xx = [ceil((totDelaySamp - sampBefore)/fs)- 1/fs : 1/fs : ceil(size(meansArr,1)/fs)] - totDelaySamp;
    xx = (1/fs : 1/fs : ceil(size(averageFluor,1)/fs)) - totDelaySamp;
    xx = xx';

    %find max dff to make a bulk fluorescense curve, center around baseline
    cdff = mean(dff,2) - blVal(fileN);
    %sem    = std(dff,0,2)./sqrt(size(dff,2));

    [maxVal, maxIdx] = max(cdff(totDelaySamp:end)); 
    peakVal(fileN) = maxVal;
    % semarr(fileN)  = sem(maxIdx);
    tmpName= strsplit(fileName,'_'); 
    if manualTitles
        concNames = TitleNames;
    else
        %concNames{ii} = tmpName{length(tmpName)-1};
        concNames{fileN} = dose;
    end

    figure
    lineProps.col{1} = [0 0.5 0];
    line([0 0],[-0.1 1])
    %mseb(xx',cdff,sem',lineProps,1);
    plot(xx,cdff)
    hold on
    try 
        plot(xx(maxIdx+totDelaySamp), maxVal,'*')
    end
    titleName = strrep(fileName,'_',' '); titleName= titleName(1:end-4);
    title(titleName)
    set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
    xlim([-totDelaySamp ceil(size(sessionImgs,3)/fs - totDelaySamp)])
    set(gca,'FontSize', 16)
    xlabel('Time (s)', 'FontSize', 22)
    ylabel('dF/F ','FontSize', 22)
    set(gca,'FontName','Arial')
    box off
    savename = fileName(1:end-4);
    %title(dose)

    %savename = [slice ligand doseOI];
    %{
    ylim([-0.1 0.3])
    xlim([-sampBefore 200])
    set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry

    set(gca,'FontSize', 16)
    xlabel('Time (s)', 'FontSize', 22)
    ylabel('dF/F ','FontSize', 22)
    %set(gcf, 'Position',  [100, 100, 600, 600])
    set(gca,'FontName','Arial')
    box off
    saveas(gcf,[savename '.svg'],'svg');
    savefig([savename '.fig']);
    %}
    if exist('varsAndFigs','dir')
        cd('varsAndFigs')
    else
        mkdir('varsAndFigs')
        cd('varsAndFigs')
    end
    saveas(gcf,[savename '.svg'],'svg');
    savefig([savename '.fig']);
    save([savename '_analyzed.mat'], 'tmpAvgFluorNoSpikes','tmpAvgFluor', 'tmpBL', 'maxVal');
    disp('Figures have been saved!')        
    cd('..')

    toc
end

    averagedArr = mean(cdff,2)';
    semArr = std(averageFluor')./sqrt(size(averageFluor,2));
    
    figure
    lineProps.col{1} = [0 0.5 0];
    line([0 0],[-0.1 1])
    if length(workFiles)>1
        mseb(xx',averagedArr,semArr,lineProps,1);
    else
        plot(xx,averagedArr)
    end
    hold on
    title(doseOI)
    xlim([-totDelaySamp ceil(size(sessionImgs,3)/fs - totDelaySamp)])

    if exist('varsAndFigs','dir')
        cd('varsAndFigs')
    else
        mkdir('varsAndFigs')
        cd('varsAndFigs')
    end
    saveas(gcf,[doseOI '.svg'],'svg');
    savefig([doseOI '.fig']);
    save([doseOI '.mat'], 'averagedArr', 'semArr','xx','-v7.3');
    disp('Figures have been saved!')        
    cd('..')
    
end


%% extract peaks from .mat files and make concentration response curves

ligand = 'DA';
doseArr = {'1nM', '3nM','10nM', '30nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'};
ExpStruct = struct();
peakArr = [];
ExpStruct.name = [ligand]; 


slice = {}; doses = {};
date = {}; doseVals = [];
ptr = 1;

for jj = 1:length(doseArr) 
    workFiles = dir(['*' doseArr{jj} '*_Analyzed.mat']);
    doseOI = doseArr{jj};
    
    for fileN = 1:length(workFiles) 
        fileName = workFiles(fileN).name;
        splitStr = regexp(fileName,'_','split');
        len  = length(splitStr);
        doseOI = splitStr(contains(splitStr,ligand)); 
        doseOI = doseOI{1}; doseOI = doseOI(length(ligand)+1:end);
        
        doseidx = strfind(doseOI,'M'); doseNumb = doseOI(1:doseidx-2);
        magnitude = doseOI(doseidx-2:end);
        doseVal = str2num(doseNumb);
        if contains(doseOI, 'nM')
            doseVal = doseVal/1000;
        end
        doseVals(ptr) = doseVal;
        doses{ptr} = doseOI;
        slice{ptr} = [splitStr{1} '_' splitStr{2}];
        datetmp = splitStr(contains(splitStr,'22')); date{ptr} = datetmp{1};
        
        load(fileName)
        peakArr = [peakArr maxVal];   
%        ExpStruct.peaks = peakArr;
%        ExpStruct.sliceN = slice;
%        save([ExpStruct.name '_' slice '.mat'], 'ExpStruct');

        ptr = ptr + 1;
    end
end
slice = slice'; date = date'; doses = doses'; peakArr = peakArr'; doseVals = doseVals';
T = table(slice, date, doses, doseVals, peakArr)
writetable(T, 'dlight_ResponseCurve2.xlsx')

