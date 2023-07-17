%% initialize wash
% This script is to analyze 2p slice sessions where a wash was done
%assumes there is a folder where delays are stored
clear all


addpath('C:\Users\bruchasadmin\Documents');
addpath(genpath('G:\code\'))
pathToFile = 'G:\2P-slice exps\Sensor Control Exps\GRABDA\216740-1\washes';
cd(pathToFile)

stimON    = 0; %this script is Wash only
avgON     = 1; rollingAvgLen = 3;

%only need to change if you want to put manual title names.
manualTitles = 0;
TitleNames = {};%{'100uM 15s_2','100uM 15s 0001','100uM 15s 0002'} ;%{'Start','End'};

%do you want to reanalyze anything you already analyzed? 1 for yes
reanalyzeQ = 1;


close all 
% Can loop through fileOI 

%dose = {'100uM_15s','100uM_15s_0001','100uM_15s_0002'} 

% TO ADD: - a way to not do all files that come up in workFiles
%         - a way to add in specific delays for each session(import data from
%         csv? (handDelay)
%       -

ligand = 'DA'; %beta
%doses = {'3nM', '10nM', '30nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'};
doses = {'100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'};
%doses = {'20hz', '5hz'};
%doses = %{'100uM'};
extra_info = '';%'420s'; % can make it blank like: {''};

downsampleQ = 0;
for dd = 1:length(doses)
    filePattern = [ligand doses{dd}];
    if isempty(extra_info)
        workFiles = dir(['*' filePattern '*.mat']);    
    else
        workFiles = dir(['*' filePattern  '*' extra_info '*.mat']);
    end

    % if you don't want to reanalyze your already analyzed data
    if ~reanalyzeQ 
        workFiles = workFiles(~endsWith({workFiles.name}, 'analyzed.mat'));
    end

    blVal     = zeros(length(workFiles),1);  
    peakVal   = zeros(length(workFiles),1);
    semarr    = peakVal;
    concNames = {}; 
    tic
    
for fileN = 1:length(workFiles)
     clearvars -except dd fileN semarr concNames peakVal blVal workFiles ...
         filePattern doses ligand pathToFile manualTitles averageFluor downsampleQ reanalyzeQ cdff ...
         xx doseOI totDelaySamp sessionImgs fs extra_info TitleNames rollingAvgLen avgON
     
    fileName = workFiles(fileN).name;

    %confirm the filename is the ligand of interest
    if isempty(strfind(fileName,ligand)) 
        continue
    end

    splitStr = regexp(fileName,'_','split');
    len  = length(splitStr);
    doseidx = find(contains(splitStr,ligand)==1); doseidx = doseidx(end);
    doseOI = splitStr(doseidx); 
    doseOI = doseOI{1}; dose = doseOI(length(ligand)+1:end);
    slice = splitStr{2};
%     date = splitStr(contains(splitStr,'22')); date = date{1};

    % load delay and file
    %CUSTOMIZE
    %load([pathToFile(1) ':\2P-slice exps\NE DA proj\delays\delay' date '.mat']); %get delay
    delay=150; handDelay = 0;
    load(fileName)
    
    if fileN == 1 % initialize
        averageFluor = zeros(size(sessionImgs,3),length(workFiles)); %CUSTOMIZE %size(sessionImgs,3),length(workFiles));
    end

    fs = 1/1.088;
    warning('do we need to correct the fs?')
    %{
    laterDate = compareDates(fileName, '220214'); % '220214' is when i switched to resonant
    if laterDate
        fs = 3; %frame avging of 10hz
        artifactThreshold = 300; 
    else
        fs = 1/1.088;
        artifactThreshold = 500; 
    end
    dsFactor = 3/1.088;
%}

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

    if avgON
        filteredImgs = imboxfilt3(sessionImgs,[1 1 rollingAvgLen]);
    end
 
    %find sum of fluorescence per area
    pixels = sum(ROImask(:));
    contourImgs = (filteredImgs.*ROImask)/pixels; 
    % find fluorescence change per area
    tmpAvgFluor = squeeze(sum(contourImgs,1)); tmpAvgFluor = squeeze(sum(tmpAvgFluor,1));
    
    if downsampleQ 
    if length(tmpAvgFluor)> 300
        sample2End = ceil(size(sessionImgs,3)/3);
        if sample2End >276
            sample2End = 276;
        tmpAvgFluor = decimate(tmpAvgFluor,3); tmpAvgFluor= tmpAvgFluor(3:sample2End+2);
        elseif sample2End < 276 && size(sessionImgs,3)>3*276
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
    end
    tmpAvgFluorNoSpikes = medfilt1(tmpAvgFluor',3);
    f0  = mean(tmpAvgFluorNoSpikes(baselineSamples), 1); f0 = mean((f0(:))~=0);
    
    % correct sizes of arrays
    if length(tmpAvgFluorNoSpikes)>length(averageFluor)
        tmpAvgFluorNoSpikes = tmpAvgFluorNoSpikes(1:length(averageFluor),:);
        averageFluor(:,fileN) = tmpAvgFluorNoSpikes;
    elseif length(tmpAvgFluorNoSpikes)<length(averageFluor)
        averageFluor = averageFluor(1:length(tmpAvgFluorNoSpikes),:);
        averageFluor(:,fileN) = tmpAvgFluorNoSpikes;
    else
        averageFluor(:,fileN) = tmpAvgFluorNoSpikes;
    end

    % Calculate df/ft
    dff = (tmpAvgFluorNoSpikes - f0)./f0;
    if size(dff,1) == 1
        dff = dff';
    end
    %[max_dff, max] = max(dff);
    tmpBL = mean(dff(baselineSamples,:));
    blVal(fileN) =tmpBL;
    % Create a time array
    %xx = [ceil((totDelaySamp - sampBefore)/fs)- 1/fs : 1/fs : ceil(size(meansArr,1)/fs)] - totDelaySamp;
    xx = (1/fs : 1/fs : ceil(size(averageFluor,1))/fs) - totDelaySamp/fs;
    xx = xx';

    %find max dff to make a bulk fluorescense curve, center around baseline
    cdff = mean(dff,2) - blVal(fileN);
    
    DurToCenter = 10; %seconds
    x0 = find(xx>0); x0tocentDur = x0 + ceil(DurToCenter*fs);
    cdff = cdff - mean(cdff(x0:x0tocentDur));
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


     %% Create heat map from before to after stim

    figure
    imagesc(1000*contourImgs(:,:,maxIdx)-contourImgs(:,:,maxIdx))
    %L = line([0 0],[0 size(allData2,1)+1]);
    %set(L,'Color','white')
    %set(L,'LineWidth',1)
    %ytick(1:size(LickTrig,1))
    titleName = strrep(fileName,'_',' '); titleName= titleName(1:end-4);
    title(titleName)
    cb = colorbar;
    set(gca,'yticklabel',{[]})
    set(gca,'xticklabel',{[]})
    title(cb,'dF/F')
    caxis([-1 2]/1000)

    %% regular avg
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
    xlim([-totDelaySamp/fs ceil(size(sessionImgs,3)/fs - totDelaySamp/fs)])
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
    save([savename '_analyzed.mat'], 'tmpAvgFluorNoSpikes','tmpAvgFluor', 'tmpBL', 'maxVal', 'maxIdx');
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
    xlim([-totDelaySamp/fs ceil(size(sessionImgs,3)/fs - totDelaySamp/fs)])

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


%% average timelocked on peak

analyzedImages = dir('*DA30*analyzed*.mat');
maxIdxs = []; timeSeries = [];
for j = 1:length(analyzedImages)
    load(analyzedImages(j).name)
    
    % maxIdx should be where we focus all analysis. add all to a vector,
    % add all time series to an array
    maxIdxs = [maxIdxs maxIdx];
    timeSeries(j,:) = tmpAvgFluorNoSpikes;
    
end

fs = 3;
delay = 150; handDelay = 0;
totDelay = delay + handDelay;
totDelaySamp = ceil(totDelay*fs);
%{
% go through, and mean pad the arrays and average them to find a curve
for j = 1:size(timeSeries,1)
    %if the max value after the wash starts is prior to the avg max peak
    subvalTMP = subval(j); abssubTMP = abssub(j);
    idxNeg = subvalTMP < 0;
    % need to consider if the pos is greatest, or the neg is greatest
    if subvalTMP == maxNeg && NegGreatest
        tsNew(j, 1:abssubTMP ) = mean(timeSeries(j,1:idxs2avg));
        tsNew(j,abssubTMP+1:end) = timeSeries(j,:);
        %{
    elseif subvalTMP ~= maxNeg && subvalTMP ~= maxPos
        halfDistance = round(abssubTMP/2); 
        remainder = mod(abssubTMP,2);
        
        tsNew(j, 1:(halfDistance+remainder)-1) = mean(timeSeries(j,1:idxs2avg));
        tsNew(j, halfDistance+remainder:end-halfDistance-1) = timeSeries(j,:);
        tsNew(j, end - halfDistance:end ) = mean(timeSeries(j,end-idxs2avg:end));        
   %}
    elseif subvalTMP ~= maxNeg && subvalTMP ~= maxPos 
        oppShift = maxShift - subvalTMP;
        if idxNeg
            tsNew(j, 1:subvalTMP) = mean(timeSeries(j,1:idxs2avg));
            tsNew(j, subvalTMP+1:length(timeSeries)+subvalTMP) = timeSeries(j,:);
            tsNew(j, length(timeSeries)+subvalTMP+1:end ) = mean(timeSeries(j,end-idxs2avg:end));        
        else % if positive
            tsNew(j, 1:oppShift) = mean(timeSeries(j,1:idxs2avg));
            tsNew(j, oppShift+1:end-subvalTMP) = timeSeries(j,:);            
%            tsNew(j, length(tsNew)-subvalTMP-length(timeSeries):end-subvalTMP-1) = timeSeries(j,:);
            tsNew(j, end - subvalTMP+1:end) = mean(timeSeries(j,end-idxs2avg:end)); 
        end
    elseif subvalTMP == maxPos && ~NegGreatest
        tsNew(j, end - abssubTMP+1:end ) = mean(timeSeries(j,end-idxs2avg:end));
        tsNew(j,1:end-abssubTMP) = timeSeries(j,:);
    elseif subvalTMP == maxPos && NegGreatest
        startShift = abs(subval(maxNegIdx)) - subval(maxPosIdx);
        tsNew(j, 1:startShift) = mean(timeSeries(j,1:idxs2avg));
        tsNew(j, startShift+1:length(timeSeries)+startShift) = timeSeries(j,:);
        tsNew(j, length(timeSeries)+startShift+1:end ) = mean(timeSeries(j,end-idxs2avg:end));        
    end
end
%}

timeSeries2 = movmean(timeSeries,15,2);

[maxval maxIdx] = max(maxIdxs); %maxIdx will not move, but the max# of idx will be padded at end
[minval minIdx] = min(maxIdxs); %minIdx will move most
maxShift = maxval-minval;
shifts  = maxIdxs - maxval;

tsNew = zeros(size(timeSeries,1),maxShift+size(timeSeries,2));
idxs2avg = 10;

blSamps = totDelaySamp-30*fs:totDelaySamp;
timeSeries = timeSeries - mean(timeSeries(blSamps));

% shifting all of the timeseries to the right
for j = 1:size(timeSeries,1)    

    %if the max value after the wash starts is prior to the avg max peak
    if j == minIdx
        % the earliest peak will shift by the value of the largest peak
        % minus the smallest peak
        tsNew(j, 1:maxShift) = mean(timeSeries(j,1:idxs2avg));
        tsNew(j, maxShift+1:maxShift+length(timeSeries)) = timeSeries(j,:);
    elseif j == maxIdx
        tsNew(j, 1:length(timeSeries)) = timeSeries(j,:); 
        tsNew(j, end-maxShift:end) = mean(timeSeries(j,end-idxs2avg:end));
    else %intermediary shifts
        initShift = abs(shifts(j));
        leftoverShift = maxShift - initShift; %for the end 
        tsNew(j, 1:initShift) = mean(timeSeries(j,1:idxs2avg));
        tsNew(j, initShift+1:length(timeSeries)+initShift) = timeSeries(j,:);        
        tsNew(j, end-leftoverShift+1:end) = mean(timeSeries(j,end-idxs2avg:end));
    end    
  
end





%% check ur work

blSamps = totDelaySamp-30*fs:totDelaySamp;
tsNew = tsNew - mean(tsNew(:,blSamps),2);

figure
for iii = 1:3
hold off
subplot(3,1,iii)
plot(timeSeries(iii,:))
hold on
plot(maxIdxs(iii)+totDelaySamp-1,timeSeries(iii,maxIdxs(iii)+totDelaySamp-1),'*')
end
sgtitle('three unshifted plots of 100uM washes')

%do we want to do a moving average of the timeseries?
timeSeries2 = movmean(timeSeries,3,2);
figure
for iii = 1:3
hold off
subplot(3,1,iii)
plot(timeSeries2(iii,:))
hold on
plot(maxIdxs(iii)+totDelaySamp-1,timeSeries2(iii,maxIdxs(iii)+totDelaySamp-1),'*')
end
sgtitle('three moving avg plots of 100uM washes')

figure
for iii = 1:3
hold off
subplot(3,1,iii)
plot(tsNew(iii,:))
hold on
plot(maxShift+minval+totDelaySamp-1,tsNew(iii,maxShift+minval+totDelaySamp-1),'*')
end
sgtitle('three plots aligned by peak')

% another way to check
[a b] = max(tsNew');

%% find average
xx = 1/fs:1/fs:length(tsNew)/fs;
tsNew = tsNew - mean(tsNew(blSamps));

sem = std(tsNew,0,1)./sqrt(size(tsNew,1)); % sem = std/sqrt(n)
lineProps.col{1} = [0 0.5 0];

figure
tsNew2 = movmean(tsNew,4,2); % MOVING AVERAGE TO SMOOTH
mseb(xx,mean(tsNew2),sem,lineProps,1);

  max(mean(tsNew))

%% extract peaks from .mat files and make concentration response curves


pwdir = pwd;
if ~contains(pwdir,"varsAndFigs")
    cd("varsAndFigs")
end

%ligand = 'Met';
ligand = 'DA'; %beta
%doses = {'3nM', '10nM', '30nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'};
doseArr =  {'30nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'};
%doseArr = {'1nM', '10nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'};
%doseArr = {'10uM'}
ExpStruct = struct();
peakArr = [];
ExpStruct.name = [ligand]; 


slice = {}; doses = {};
date = {}; doseVals = [];
ptr = 1;

for jj = 1:length(doseArr) 
    workFiles = dir(['*' ligand doseArr{jj}  '*_analyzed.mat']);
    doseOI = doseArr{jj};
    
    for fileN = 1:length(workFiles) 
        fileName = workFiles(fileN).name;
        splitStr = regexp(fileName,'_','split');
        len  = length(splitStr);
        %doseOI = splitStr(contains(splitStr,ligand)); 
        %doseOI = doseOI{1}; doseOI = doseOI(length(ligand)+1:end);
        
        %doseidx = strfind(doseOI,'M'); doseNumb = doseOI(1:doseidx-2);
        %magnitude = doseOI(doseidx-2:end);
        doseNumb = doseOI(1:end-2);
        doseVal = str2num(doseNumb);
        if contains(doseOI, 'nM')
            doseVal = doseVal/1000;
        end
        doseVals(ptr) = doseVal;
        doses{ptr} = doseOI;
        slice{ptr} = [splitStr{1} '_' splitStr{2}];
        %datetmp = splitStr(contains(splitStr,'22')); date{ptr} = datetmp{1};
        
        load(fileName)
        peakArr = [peakArr maxVal];   
%        ExpStruct.peaks = peakArr;
%        ExpStruct.sliceN = slice;
%        save([ExpStruct.name '_' slice '.mat'], 'ExpStruct');

        ptr = ptr + 1;
    end
end
slice = slice'; %date = date'; 
doses = doses'; peakArr = peakArr'; doseVals = doseVals';
%T = table(slice, date, doses, doseVals, peakArr)
T = table(slice, doses, doseVals, peakArr)
writetable(T, [ligand '_ResponseCurve2.xlsx'])

