% things to add
% automatic detection of when stim is on  
clear all
close all
clc

codeDir = 'C:\Users\avima\OneDrive\Documents\MATLAB\Starter Code - photometry'; %CUSTOMIZE
dataDir = 'G:\2P-slice exps\Sensor Control Exps\GRABDA\216740-1\stim\saved ROIs'; %CUSTOMIZE
addpath(codeDir);
addpath(genpath('G:\code\'));

cd(dataDir)

stimLength = [3 30]; %[3 15 30] in seconds

%fs for galvano scan
fs = 1/1.088; % scans per second
details.rollingAvgLen = 3;
delay = 15; %in seconds

% if rolling averaged:
%12
%19 44

details.preStimSamp = 13; %floor(delay*fs)-1;
details.postStimSamp= [19 44]; % [19 30 44] ceil(delay*fs+rollingAvgLen-1) + ceil(stimLength*fs);

details.manualTitles = 0;
details.saveON = 0;
details.avgON = 1;


details.downsampleQ = 0;
details.freq = 5; details.stimTime = 3; sessionToSave = [ num2str(details.freq) 'hz_' num2str(details.stimTime) 's'];
workFiles = dir(['*' num2str(details.freq) 'hz_' num2str(details.stimTime) 's*']);
%{
%dirOI = uigetfile('*.csv');
slice = 'S2';
ligand = 'Leu';
dose = { '100uM'};%{'1nM', '10nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'} 

%only need to chnage if you want to put manual title names.
manualTitles = 1;
TitleNames = {'Start','End'};
%}
%% create data array and avg


subjects = {}; subPtr = 1; trials = [];
for fileN = 1:length(workFiles)
     clearvars -except dd fileN semarr concNames peakVal blVal workFiles delay details...
         filePattern doses ligand pathToFile averageFluor reanalyzeQ cdff slices ...
         xx doseOI totDelaySamp sessionImgs fs currTrials subjects subPtr trials trialsCell subjectsCell slicesCell stimTime
     
    fileName = workFiles(fileN).name;
    splitStr = regexp(fileName,'_','split');
    currSubject  = splitStr{1};
    currSlice    = splitStr{2};
    if ~isempty(subjects)
        if ~strcmp(subjects{subPtr},currSubject)
            subPtr = subPtr + 1;
            subjects{subPtr} = currSubject;
            trials{subPtr} = currTrials;
            slices{subPtr} = currSlice;
            currTrials = 1; 
        else %Subjects exists, and we are not on the first trial
            currTrials = currTrials + 1; %add one 
            trials{subPtr} = currTrials;
            slices{subPtr} = currSlice;            
        end
    else % first subject
        subjects = {currSubject};
        currTrials = 1;
    end

    %has each trial and cell written out each time 
    subjectsCell{fileN} = subjects{subPtr};
    trialsCell{fileN} = currTrials;
    slicesCell{fileN} = currSlice;


    len  = length(splitStr);
    delayOI = splitStr(contains(splitStr,'delay')); 
    slice = splitStr{2};
    load(fileName)
    
    if fileN == 1 % initialize
        averageFluor = zeros(size(sessionImgs,3),length(workFiles)); %CUSTOMIZE %size(sessionImgs,3),length(workFiles));
    end
    
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
    totDelay = delay; totDelaySamp = ceil(totDelay*fs); %the total delay in samples

    %samples before for baseline
    sampBefore = ceil(15*fs); %30 s
    baselineSamples = totDelaySamp - sampBefore+1:totDelaySamp; 
 
    if exist('sessionImgs','var')
        meanBL = mean(sessionImgs(:,:,baselineSamples),3); meanBL = imgaussfilt(meanBL,2);
    else
        warning([fileName ' has no sessionImgs saved!!!'])
        continue
    end
    
    if details.avgON
        filteredImgs = imboxfilt3(sessionImgs,[1 1 details.rollingAvgLen]);
    end

    %find sum of fluorescence per area
    pixels = sum(ROImask(:));
    contourImgs = (filteredImgs.*ROImask)/pixels; 
    
    % find fluorescence change per area
    tmpAvgFluor = squeeze(sum(contourImgs,1)); tmpAvgFluor = squeeze(sum(tmpAvgFluor,1));
    
    sesSamps    = length(tmpAvgFluor);  
    sample2End = sesSamps; %ceil(size(sessionImgs,3)/3);

    if details.downsampleQ 
        if length(tmpAvgFluor)> 120 %120 bc size of 2 min session in galvo is 114
                tmpAvgFluor = decimate(tmpAvgFluor, 3); 
                tmpAvgFluor = tmpAvgFluor(1:sample2End); 
        end
    end
    %{
if length(tmpAvgFluor)> 120 %120 bc size of 2 min session in galvo is 114
        sample2End = size(averageFluor,1); %ceil(size(sessionImgs,3)/3);
        if sample2End > sesSamps
            sample2End = sesSamps;
        tmpAvgFluor = decimate(tmpAvgFluor,3); tmpAvgFluor= tmpAvgFluor(3:sample2End);
        elseif sample2End < sesSamps 
             tmpAvgFluor = decimate(tmpAvgFluor, 3); 
             bltmp(1:sample2End) =  mean(tmpAvgFluor(1:sample2End),1);
             if size(tmpAvgFluor,1)==1
                 tmpAvgFluor = tmpAvgFluor';
             end
            tmpAvgFluor= [bltmp'; tmpAvgFluor(1:sample2End)];
        end
    elseif size(sessionImgs,3)<sesSamps
         if size(tmpAvgFluor,1)==1
                 tmpAvgFluor = tmpAvgFluor';
         end
        bltmp(1:sesSamps-length(tmpAvgFluor)) =  mean(tmpAvgFluor(1:sesSamps-length(tmpAvgFluor)),1);
        tmpAvgFluor= [bltmp'; tmpAvgFluor];            
    end
    end
    %}
    
    
    tmpAvgFluorNoSpikes = medfilt1(tmpAvgFluor',3);
    f0  = mean(tmpAvgFluorNoSpikes(baselineSamples), 1); f0 = mean((f0(:))~=0);
    
    
    %averageFluor(:,fileN) = tmpAvgFluorNoSpikes;
    
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
    xx = (1/fs : 1/fs : length(tmpAvgFluor)/fs) - totDelaySamp/fs;
    xx = xx';

    %find max dff to make a bulk fluorescense curve, center around baseline
    cdff = mean(dff,2) - blVal(fileN);
    %sem    = std(dff,0,2)./sqrt(size(dff,2));

    [maxVal, maxIdx] = max(cdff(totDelaySamp:end)); 

    peakVal(fileN) = maxVal;
    % semarr(fileN)  = sem(maxIdx);
    tmpName= strsplit(fileName,'_'); 
    
    if details.stimTime == 3
        if sesSamps == 180
            preSamps = 37:45;
            postSamps= 55:63;
        elseif sesSamps < 120
            preSamps = 10:13;
            postSamps= 18:21;
        end
    else
        if fs == 3
            preSamps = 37:45;
            postSamps= 55:63;
        elseif fs == 1/1.088
            preSamps = 10:13;
            postSamps= 43:46;
        end
    end
    
    prePostVal = [mean(cdff(preSamps)) mean(cdff(postSamps))];


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
    
    if details.saveON
        saveas(gcf,[savename '.svg'],'svg');
        savefig([savename '.fig']);
        save([savename '_analyzed.mat'], 'tmpAvgFluorNoSpikes','tmpAvgFluor', 'tmpBL', 'prePostVal', 'currSubject', 'trials','subjects','fs');
        disp('Figures have been saved!')        
    end
    cd('..')

end

%% Save pre and post means
peakVals = [];
preVals = [];
cd('varsAndFigs')

meansDir = dir(['*' num2str(details.freq) 'hz_' num2str(details.stimTime) 's_*_analyzed.mat']);
sessionToSave = [ num2str(details.freq) 'hz_' num2str(details.stimTime) 's'];
for i= 1:length(meansDir)
    load(meansDir(i).name)    
    % find relevant info on trials
    trialNum{i} = num2str(trialsCell{i});
    peakVals = [peakVals prePostVal(2)];
    preVals = [preVals prePostVal(1)];
        
end

subNames2 = char(subjectsCell);
slices2   = char(slicesCell);
trials2   = char(trialNum);
peakVals  = reshape(peakVals,[],length(subjectsCell))';
preVals = reshape(preVals,[],length(subjectsCell))';
T = table(subNames2, slices2, trials2, peakVals,preVals);
newsavename = [sessionToSave '_peaks.csv'];
writetable(T, newsavename)
cd('..')

%% AUC + Peak analysis ?
cd('varsAndFigs')

meansDir = dir('*_analyzed.mat');


for j = 1:length(meansDir)
    load(meansDir(j).name)
% save variables
trials   = trials;
subNames = subjects;
titleName = titleName;
savename = savename;
%fs = dataForPlots.FS;
xlims = [-15 45];
timeBefore = totDelaySamp;

% make uniform sizing across data arrays and time vector
if exist('timevec')
    if size(timevec,1) == 1
    timevec = timevec';
    end
else
    if exist('time')
        timevec = time;
        if size(timevec,1) == 1
            timevec = timevec';
        end
    else
        warning('you might not have a time vector')
        timevec = xx;
    end
end

% define eventEnd (in s) for the sake of peak extraction and AUC 
eventEnd = 30;
endSamp = round((eventEnd + timeBefore)*fs);
startSamp = round(timeBefore*fs);

cTrials  = cumsum(trials);
dcTrials = cTrials - cTrials(1)+1; % start of each trial

names = []; 
peakVals = []; AUCarr = [];
preVals = [];
for i = 1:length(trials)
    
    % find relevant info on trials
    trialAmt = trials(i); 
    trialsOI = dcTrials(i):cTrials(i);
    currTime = timevec(startSamp:endSamp);
    currData = allData(startSamp:endSamp,trialsOI);
    currTrial = find(trialsOI == i);
    
    AUC = trapz(currTime, currData);    
    [maxVal,maxIdx] = max(currData);

    peakVals = [peakVals prePostVal(2)];
    preVals = [preVals prePostVal(1)];
    AUCarr = [AUCarr AUC];
        
end
end
subNames2 = char(subNames');
peakVals  = reshape(peakVals,[],length(subNames))';
preVals = reshape(preVals,[],length(subNames))';
AUCarr    = reshape(AUCarr,[],length(subNames))';
T = table(subNames2, peakVals,preVals, AUCarr);
newsavename = [savename '_AUC.csv'];
writetable(T, newsavename)



%% Avg trace?





%% OLD

%%

%allDelays = [20 20 20 20 20 20 20 20];% [3 4 4 6 4 8 3 7]; %delay from starting session to switching tubes to new agent
%defaultDelay = 105; %it takes 1 min, 45 seconds for the pump to send agents
stimStruct = struct(); 
stimValues = [];

% Can loop through fileOI 
for j = 1:length(stimLength)
    stimOI = stimLength(j);    
    titleNames{j} = ['Chrimson stim for' num2str(stimLength(j)) ' s'];
    dirOI = dir(['*chrim*' slice '*' num2str(stimOI) 's_*.csv']);
    savename = ['chrim_' slice '_' num2str(stimOI) 's'];
    peaks = []; concNames = {}; semarr = []; preStimVal = [];
    postStimVal= [];
for ii = 1:length(dirOI)
    fileOI = dirOI(ii).name;
    t = readtable(fileOI);

%find which columns are the means of ROI
colNames = t.Properties.VariableNames;
meanIdx  = contains(colNames,'Mean');
meansTable    = t(:,meanIdx);

% Calculate df/f
meansArr = table2array(meansTable);
f0  = mean(meansArr(details.preStimSamp)); %mean(meansArr(1:preStimSamp));
%dff = (meansArr(totDelaySamp - sampBefore:end,:) - f0)./f0;
dff = (meansArr - f0)./f0;
%[max_dff, max] = max(dff);

% Create a time array
%xx = [ceil((totDelaySamp - sampBefore)/fs)- 1/fs : 1/fs : ceil(size(meansArr,1)/fs)] - totDelaySamp;
xx = [1/fs : 1/fs : ceil(size(meansArr,1)/fs)]-details.preStimSamp;
xx = xx';

%find max dff to make a bulk fluorescense curve
avgdff = mean(dff,2);
sem    = std(dff,0,2)./sqrt(size(dff,2));

%save preStim value and postStim value
preStimVal = [preStimVal dff(details.preStimSamp)]; postStimVal = [postStimVal dff(details.postStimSamp(j))]; 

%save peaks
%
[maxVal, maxIdx] = max([avgdff(1:details.preStimSamp); avgdff(details.postStimSamp(j):end)]); %adjusted 120621
if maxIdx>details.postStimSamp
    maxIdx = maxIdx + ceil(delay*fs+details.rollingAvgLen-1) + ceil(stimLength(j)*fs);
end
peaks = [peaks maxVal];
%semarr= [semarr sem(maxIdx)];
tmpName= strsplit(fileOI,'_'); 

%figure
%lineProps.col{1} = [0 0.5 0];
%mseb(xx',avgdff,sem',lineProps,1);
%line([0 0],[-1 1])
%line([postStimSamp(j)-preStimSamp postStimSamp(j)-preStimSamp],[-1 1])
%plot(xx,avgdff)
%hold on
try 
%    plot(xx(maxIdx+totDelaySamp), maxVal,'*')
end
%title(titleNames{j})
end
save([savename '.mat'], 'preStimVal', 'postStimVal');

end
%% 

stimStruct = setfield(stimStruct, ['chrim_' num2str(stimOI) 's'], stimValues);

savename = [slice stimOI];

ylim([-0.5 0.5])
xlim([-details.preStimSamp 200])
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry

set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('dF/F ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off
saveas(gcf,[savename '.svg'],'svg');
savefig([savename '.fig']);
save([savename '.mat'], 'peaks', 'semarr');
%}

%% extract peaks from .mat files and make concentration response curves


slice = 'S2';
ligand = 'beta';
dose = {'1nM', '10nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'} 
ExpStruct = struct();
peakArr = [];
ExpStruct.name = [ligand]; 

for jj = 1:length(dose) 
    stimOI = dose{jj};
    load([slice ExpStruct.name stimOI '.mat'])
    peakArr = [peakArr peaks];   
    ExpStruct.peaks = peakArr;
    
end