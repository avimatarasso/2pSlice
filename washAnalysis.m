%% initialize wash
% This script is to analyze 2p slice sessions where a wash was done
%assumes there is a folder where delays are stored
clear all

addpath('F:\code\2pSliceAnalysis\oir2stdData-master');
addpath('C:\Users\avima\OneDrive\Documents\GitHub\2pSlice');
pathToFile = 'F:\2P-slice exps\NE DA proj\GRABNE\washes\';%'F:\2P-slice exps\022222\';%'F:\2P-slice exps\022222\105304-3_S1_dL_CA1_20hz_30s_15sdelay.oir';
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

ligand = 'NE'; %beta
doses = {'3nM', '10nM', '30nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'};
for dd = 1:5%:length(doses)
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
    
    load(['F:\2P-slice exps\NE DA proj\delays\delay' date '.mat']); %get delay
    load(fileName)
    if fileN == 1 % initialize
        averageFluor = zeros(size(CRGB,4),length(workFiles));
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

    % totDelaySamp = sampBefore+1; %%CHANGE
   % averageFluor = zeros(size(CRGB,4),length(workFiles));
    %semArr = zeros(size(CRGB,4),1);

    baselineSamples = totDelaySamp - sampBefore:totDelaySamp; 
    if exist('CRGB','var')
        meanBL = mean(CRGB(:,:,:,baselineSamples),4); meanBL = imgaussfilt(meanBL,2);
    else
        warning([fileName ' has no CRGB saved!!!'])
        continue
    end
    a_gray = rgb2gray(meanBL); level = 0.1;
    a_gray = imbinarize(a_gray,level); 
    pixels = sum(a_gray(:));

    % a_bw is a thresholded contour of your region. 
    % ADD: sanity check it before calculating avg
    tmpImg = zeros(size(CRGB,1),size(CRGB,2));
   contourImgs = zeros(size(CRGB,1),size(CRGB,2), size(CRGB,4));
   tmpAvgFluor = zeros(size(CRGB,4),1);
    for imgNumb = 1:size(CRGB,4) 
        tmpImg = rgb2gray(CRGB(:,:,:,imgNumb));
        tmpImg = immultiply(tmpImg,a_gray); contourImgs(:,:,imgNumb) = tmpImg;
        tmpAvgFluor(imgNumb) = sum(tmpImg(:))/pixels;
        %semArr(imgNumb) = std(tmpImg(:))./sqrt(size(tmpImg,2)); % ADD: do i want to divide by area?? 
    end
    f0  = mean(contourImgs(:,:,baselineSamples),3); f0 = mean((f0(:))~=0);
    if length(tmpAvgFluor)> 300
        sample2End = ceil(size(CRGB,4)/3);
        if sample2End >276
            sample2End = 276;
        tmpAvgFluor = decimate(tmpAvgFluor,3); tmpAvgFluor= tmpAvgFluor(3:sample2End+2);
        elseif sample2End < 276 && size(CRGB,4)>276
             tmpAvgFluor = decimate(tmpAvgFluor,3); 
             bltmp(1:276 - sample2End) =  mean(tmpAvgFluor(1:276 - sample2End),1);
            tmpAvgFluor= [bltmp'; tmpAvgFluor(1:sample2End)];
        end
    elseif size(CRGB,4)<276
         
        bltmp(1:276-length(tmpAvgFluor)) =  mean(tmpAvgFluor(1:20),1);
        tmpAvgFluor= [bltmp'; tmpAvgFluor];            
    end
    tmpAvgFluorNoSpikes = medfilt1(tmpAvgFluor,3);
    averageFluor(:,fileN) = tmpAvgFluorNoSpikes;
    % adjust averageFluor to be the same size decimate
    %
    %{
    if fileN > 1 && length(find(tmpAvgFluor~=0)) < length(find(tmpAvgFluor~=0)) %CUSTOMIZE
        tmpBef = averageFluor(:,fileN-1); tmp2Idx = find(tmpBef==0); tmpBef(tmp2Idx) = [];  
        lenTmpBef = length(tmpBef);
        tmpIdx = find(tmpAvgFluor==0); tmpAvgFluor(tmpIdx) = [];  
        lenTmp = length(tmpAvgFluor);
        
        tmpBef = decimate(tmpBef,3); tmpBef= tmpBef(3:lenTmp+2);
    elseif fileN > 1 && length(find(tmpAvgFluor~=0)) > length(find(tmpAvgFluor~=0))
        tmpBef = averageFluor(:,fileN-1); tmp2Idx = find(tmpBef==0); tmpBef(tmp2Idx) = [];  
        lenTmpBef = length(tmpBef);
        tmpIdx = find(tmpAvgFluor==0); tmpAvgFluor(tmpIdx) = [];  
        lenTmp = length(tmpAvgFluor);
        
        tmpAvgFluor = decimate(tmpAvgFluor,3); tmpAvgFluor= tmpAvgFluor(3:lenTmpBef+2);
    end
    if fileN > 3
            averageFluor = [averageFluor(:,1:fileN-2) tmpBef tmpAvgFluor];
    elseif fileN == 2
        averageFluor = [tmpBef tmpAvgFluor];
    end
    averageFluor(imgNumb,fileN)
%}
    
    % Calculate df/f
    dff = (averageFluor(:,fileN) - f0)./f0;
    %[max_dff, max] = max(dff);
    tmpBL = mean(dff(baselineSamples,:));
    blVal(fileN) =tmpBL;
    % Create a time array
    %xx = [ceil((totDelaySamp - sampBefore)/fs)- 1/fs : 1/fs : ceil(size(meansArr,1)/fs)] - totDelaySamp;
    xx = [1/fs : 1/fs : ceil(size(averageFluor,1)/fs)]-totDelaySamp;
    xx = xx';

    %find max dff to make a bulk fluorescense curve, center around baseline
    avgdff = mean(dff,2) - blVal(fileN);
    %sem    = std(dff,0,2)./sqrt(size(dff,2));

    [maxVal, maxIdx] = max(avgdff(totDelaySamp:end)); %adjusted 120621
    peakVal(fileN) = maxVal;
%    semarr(fileN)  = sem(maxIdx);
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
    %mseb(xx',avgdff,sem',lineProps,1);
    plot(xx,avgdff)
    hold on
    try 
        plot(xx(maxIdx+totDelaySamp), maxVal,'*')
    end
    titleName = strrep(fileName,'_',' '); titleName= titleName(1:end-4);
    title(titleName)
    set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
    xlim([-totDelaySamp ceil(size(CRGB,4)/fs - totDelaySamp)])
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
    save([savename '_Analyzed.mat'], 'tmpAvgFluorNoSpikes','tmpAvgFluor', 'tmpBL', 'maxVal');
    disp('Figures have been saved!')        

    cd('..')

    toc
end
    
    averagedArr = mean(averageFluor,2)';
    semArr = std(averageFluor')./sqrt(size(averageFluor,2)); 
    
    figure
    lineProps.col{1} = [0 0.5 0];
    line([0 0],[-0.1 1])
    mseb(xx',avgdff,semArr,lineProps,1);
    %plot(xx,averagedArr)
    hold on
    title(doseOI)
    xlim([-totDelaySamp ceil(size(CRGB,4)/fs - totDelaySamp)])

    if exist('varsAndFigs','dir')
        cd('varsAndFigs')
    else
        mkdir('varsAndFigs')
        cd('varsAndFigs')
    end
    saveas(gcf,[doseOI '.svg'],'svg');
    savefig([doseOI '.fig']);
    save([doseOI '_Analyzed.mat'], 'averagedArr', 'semArr','xx');
    disp('Figures have been saved!')        

    cd('..')

end


%% extract peaks from .mat files and make concentration response curves


slice = 'S2';
ligand = 'beta';
dose = {'1nM', '10nM', '100nM', '300nM', '1uM', '3uM', '10uM', '30uM', '100uM'} 
ExpStruct = struct();
peakArr = [];
ExpStruct.name = [ligand]; 

for jj = 1:length(dose) 
    doseOI = dose{jj};
    load([slice ExpStruct.name doseOI '.mat'])
    peakArr = [peakArr peaks];   
    ExpStruct.peaks = peakArr;
    ExpStruct.sliceN = slice;
    save([ExpStruct.name '_' slice '.mat'], 'ExpStruct');

end




