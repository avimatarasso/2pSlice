% This script is to analyze 2p slice sessions
%{
addpath('G:\code\2pSliceAnalysis\oir2stdData-master');
pathToFile = 'G:\2P-slice exps\022222\105304-3_S1_dL_CA1_5hz_3s_15sdelay.oir';

% extract the metadata
strfind(pathToFile,'s_');
splitStr = regexp(pathToFile,'_','split');
len = length(splitStr);
freq = splitStr{len-2};
stimTime = splitStr{len-1};

if strcmp(stimTime, '3s')
    stimLims = [46 54];
elseif strcmp(stimTime, '30s')
    stimLims = [46 135]; %need to check this
end

tic
[a,stdData]=oir2stdData(pathToFile);
b = (stdData(1).Image{1});
c = im2uint8(squeeze(b)); 

%uint16 to RGB
CRGB = zeros(size(c,1),size(c,2),3,size(c,3));
%
for imgNumb = 1:size(c,3)
    C = mat2gray(c(:,:,imgNumb));
    tmp = cat(3, C, C, C);
    CRGB(:,:,:,imgNumb) = tmp;
end
%
baselineIMG = 1:15; meanBL = mean(CRGB(:,:,:,baselineIMG),4);

meanBLb = imgaussfilt(meanBL,2);
a_gray = rgb2gray(meanBLb);
level = 0.1;
a_bw = imbinarize(a_gray,level); 
%imshow(c(:,:,15)); hold on 

% a_bw is a thresholded contour of your region. we should check it before
% calculating avg
tC = zeros(size(CRGB,1),size(CRGB,2),size(CRGB,4));
%tCArt = tC;

for imgNumb = 1:size(CRGB,4)
    tmpImg = rgb2gray(CRGB(:,:,:,imgNumb));
    tCtmp = immultiply(tmpImg,a_bw);
    tmpArt = immultiply(tmpImg,~a_bw); % for finding artifacts
    
    %IF RESONANT
    if imgNumb >= stimLims(1) && imgNumb <= stimLims(2)
    %if imgNumb>45 %starting img %TOO MUCH HARDCODDE
    %if dS > 50 %if this image has artifact ....? (Is this incorrect?)
    %if length(dS) == 1 %redundant?
    
    % find sum of rows in artifactual image - avg image
    %tmpArt = imgaussfilt(tmpArt,2);
    tmpArtS = sum(tmpArt - a_gray,2);

    [m,i]=maxk(tmpArtS,100); %Artifacts should be ~80 rows
    diffSort = find(diff(sort(i))>10); %has the diff between top rows 
    dS = diff(diffSort);

    % BAD: can i just get rid of > 3.5
    %rows2Nan = find(tmpArtS > 3);
    %for rr = 1:length(rows2Nan)
    %    tCtmp(rows2Nan(rr),:) = ones(1,size(tCtmp,2)); 
    %end
    
    %{
    %max approach
    [m,i]=maxk(tmpArtS,100); %Artifacts should be ~80 rows
    diffSort = find(diff(sort(i))>1); %has the diff between top rows 
    dS = diff(diffSort);

    %
    %sum approach
    i= (find(tmpArtS>30)); %Artifacts should be ~80 rows
    diffSort = find(diff(sort(i))>1); %has the diff between top rows 
    dS = diff(diffSort);
    %
    %find the number of artifactual rows
    %artRows = find(abs(dS - 80) < 3); dS= dS(artRows);
    %artRows = find(dS > 50); dS = dS(artRows);
    %80 is used here for 10 frames avged during resonant
       %} 
    %the amount of artifact rows should be within 3 from 80 rows
        ArtIdx = find(dS);
        ArtSt = diffSort(ArtIdx) + 1;
        ArtEn = ArtSt + dS;
        tCtmp(ArtSt:ArtEn,:) = nan; 
        tC(:,:,imgNumb) = tCtmp;
    %    
    end
    %else
    %    imgNumb
    %end
    %end
    
    tC(:,:,imgNumb) = tCtmp;
end
toc
%}
%% BETTER
% This script is to analyze 2p slice sessions
clear all

addpath('F:\code\2pSliceAnalysis\oir2stdData-master');
pathToFile = 'F:\2P-slice exps\022222\';%'F:\2P-slice exps\022222\105304-3_S1_dL_CA1_20hz_30s_15sdelay.oir';

filePattern = '105304-3_S1_dL_CA1_20hz_3s_15sdelay';

resonant = 1;
if resonant == 1
    fs = 3; %frame avging of 10hz
else
    fs = 1/1.088;
end

workFiles = dir([pathToFile '\' filePattern '*.oir']);

preVal  = zeros(length(workFiles),1);  
postVal = zeros(length(workFiles),1);

for fileN = 1:length(workFiles)
% extract the metadata

fileName = workFiles(fileN).name; 
strfind(fileName,'s_');
splitStr = regexp(fileName,'_','split');
len  = length(splitStr);
freq = splitStr(contains(splitStr,'hz')); freq= freq{1};

delayStr = splitStr(contains(splitStr,'sdelay')); delayStr = delayStr{1};
stimTimeStr = splitStr(find(contains(splitStr,'sdelay'))-1); stimTimeStr = stimTimeStr{1};
stimTime = stimTimeStr(1:findstr(stimTimeStr,'s')-1);
stimTime = str2double(stimTime);

delay = delayStr(1:findstr(delayStr,'sdelay')-1); %make it in seconds!
delay = str2double(delay);

% Preprocessing - get rid of artifacts

%find the artifactual images
    stimLims = [floor(delay*fs+1) ceil((delay+stimTime)*fs)];
    artifactThreshold = 300; 

tic

%Extract data
[a,stdData]=oir2stdData(fileName);
b = (stdData(1).Image{1});
c = im2uint8(squeeze(b)); %turn into uint8 image

% convert uint8 to RGB to visualize
CRGB = zeros(size(c,1),size(c,2),3,size(c,3));
for imgNumb = 1:size(c,3)
    C = mat2gray(c(:,:,imgNumb));
    tmp = cat(3, C, C, C);
    CRGB(:,:,:,imgNumb) = tmp;
end

if fileN == 1
    % initialize both avg and sem
    averageArr = zeros(size(CRGB,4),length(workFiles));
    semArr = zeros(size(CRGB,4),1);
end

baselineIMG = 1:stimLims(1) - 1; meanBL = mean(CRGB(:,:,:,baselineIMG),4);
meanBLb = imgaussfilt(meanBL,2);
a_gray = rgb2gray(meanBLb); level = 0.1;
a_bw_o = imbinarize(a_gray,level); 
%imshow(c(:,:,15)); hold on 

% a_bw is a thresholded contour of your region. 
% ADD: sanity check it before calculating avg
tC = zeros(size(CRGB,1),size(CRGB,2),size(CRGB,4));
tmpArtBL = zeros(size(CRGB,1),size(CRGB,2));
% Create a time array
%xx = [ceil((totDelaySamp - sampBefore)/fs)- 1/fs : 1/fs : ceil(size(meansArr,1)/fs)] - totDelaySamp;
xx = [1/fs : 1/fs : ceil(size(CRGB,4)/fs)]' - stimLims(1)/fs;

for imgNumb = 1:size(CRGB,4)
    a_bw = a_bw_o; % reinitialize a_bw
    tmpImg = rgb2gray(CRGB(:,:,:,imgNumb));
    tCtmp = immultiply(tmpImg,a_bw);
    tmpArt = immultiply(tmpImg,~a_bw); % for finding artifacts
    
    % find background
    if imgNumb < stimLims(1)
        tmpArtBL = (tmpArtBL + tmpArt);
    elseif imgNumb == stimLims(1) %when your BL is all added up, find the avg 
        tmpArtBL = tmpArtBL/stimLims(1);
        %f0  = mean(averageArr(1:stimLims(1))); % average of the average/array for the baseline
        f0  = mean(averageArr(1:stimLims(1))); % average of the average/array for the baseline
    end
    
    % IF RESONANT - have not checked galvano
    if imgNumb >= stimLims(1) && imgNumb <= stimLims(2)
        % find sum of rows in artifactual image - avg image
        tmpArt = imgaussfilt(tmpArt,2);     tmpArtBL = imgaussfilt(tmpArtBL,2); 
        tmpArt = tmpArt - tmpArtBL; %subtract baseline
        tmpArtS = sum(tmpArt,2);

        %find anywhere above threshold
        logAboveThresh = tmpArtS > artifactThreshold*mean(tmpArtBL)'; % 300 percent  baseline?
        i = find(logAboveThresh);

        %[m,i]=maxk(tmpArtS,artifactThreshold); %Artifacts should be ~80 rows
        %dS = diff(sort(i)); %difference b/w where artifacts are
        %diffSort = find(dS>50); %has the diff between top rows 

        %dS = diff(diffSort);

        %the amount of artifact rows should be within 3 from 80 rows
        %ArtIdx = find(dS);
        %ArtSt = diffSort(ArtIdx) + 1;
        %ArtEn = ArtSt + dS;
        tCtmp(i,:) = nan; %repmat(mean(tmpArtBL(i,:),2), [1, 512]); %nan 
        a_bw(i,:) = 0; %adjust area for now
       % unnecessary?: tC(:,:,imgNumb) = tCtmp;
    end
    
    tC(:,:,imgNumb) = tCtmp;
    tCtmp(isnan(tCtmp)) = []; tCtmp(tCtmp==0) = [];
    pixels = sum(a_bw(:));
    averageArr(imgNumb,fileN) = sum(tCtmp)/pixels;
    semArr(imgNumb) = std(tCtmp,0,2)./sqrt(size(tCtmp,2)); % ADD: do i want to divide by area?? 
    
    %check tC before clearing
    %implay(tC)
end

dffTmp  = (averageArr(:,fileN) - f0)./f0;

preVal(fileN)  = dffTmp(stimLims(1)-1);
postVal(fileN) = dffTmp(stimLims(2)+1);

dff(:,fileN) = dffTmp;
dffSem = std(dffTmp')./sqrt(size(dffTmp,2));%std(dff,0,2)./sqrt(size(dff,2));

figure
lineProps.col{1} = [0 0.5 0];
%dff(stimLims) = nan;
plot(xx,dffTmp);
underLocs = findstr(fileName,'_');
tmpFileName = fileName; tmpFileName(underLocs) = ' ';
title(tmpFileName)
end
toc

dffAvg = mean(dff,2);
dffSem = std(dff')./sqrt(size(dff,2));%std(dff,0,2)./sqrt(size(dff,2));
zScore = (dffAvg - mean(dffAvg(1:stimLims(1)-1)))/std(dffAvg);
%[max_dff, max] = max(dff);

%[maxVal, maxIdx] = max([dff(1:stimLims(1)); dff(stimLims(2):end)]); 
%if maxIdx>stimLims(2)
%    maxIdx = maxIdx + ceil(delay*fs+rollingAvgLen-1) + ceil(stimLength(j)*fs);
%end
%peaks = [peaks maxVal];

figure
lineProps.col{1} = [0 0.5 0];
%dff(stimLims) = nan;
mseb(xx',zScore',dffSem,lineProps,1);
title('averaged z-score trace')
%mseb(xx',dffAvg',dffSem,lineProps,1);

save(filePattern, dffAvg,dffSem, xx, preVal, postVal)

%% Look at the images in Matlab
implay(tC)
%%
implay(CRGB)

%% find avg

% Calculate df/f
meansArr = table2array(meansTable);
f0  = mean(meansArr(preStimSamp)); %mean(meansArr(1:preStimSamp));
%dff = (meansArr(totDelaySamp - sampBefore:end,:) - f0)./f0;
dff = (meansArr - f0)./f0;
%[max_dff, max] = max(dff);

% Create a time array
%xx = [ceil((totDelaySamp - sampBefore)/fs)- 1/fs : 1/fs : ceil(size(meansArr,1)/fs)] - totDelaySamp;
xx = [1/fs : 1/fs : ceil(size(meansArr,1)/fs)]-preStimSamp;
xx = xx';

%find max dff to make a bulk fluorescense curve
avgdff = mean(dff,2);
sem    = std(dff,0,2)./sqrt(size(dff,2));

%save preStim value and postStim value
preStimVal = [preStimVal dff(preStimSamp)]; postStimVal = [postStimVal dff(postStimSamp(j))]; 

%save peaks
%
[maxVal, maxIdx] = max([avgdff(1:preStimSamp); avgdff(postStimSamp(j):end)]); %adjusted 120621
if maxIdx>postStimSamp
    maxIdx = maxIdx + ceil(delay*fs+rollingAvgLen-1) + ceil(stimLength(j)*fs);
end
peaks = [peaks maxVal];
%semarr= [semarr sem(maxIdx)];
tmpName= strsplit(fileOI,'_'); 
















%%
%{
baselineIMG = 1:15; meanBL = mean(c(:,:,baselineIMG),3);%mean(CRGB(:,:,baselineIMG),4);

meanBLb = imgaussfilt(meanBL,2);
a_gray = rgb2gray(meanBLb); level = 0.1;
a_bw = imbinarize(a_gray,level); 
imshow(meanBLb); hold on

% a_bw is a thresholded contour of your region. we should check it before
% calculating avg
tCRGB = zeros(size(CRGB));
for imgNumb = 1:size(CRGB,4)
    tCRGB = immultiply(CRGB(:,:,:,imgNumb),a_bw);
end

%{
[M,c] = contourf(a_bw);

%find the contour pts that are largest
contourPts = M(2,find(M(1,:) == level));

%finding the max k number of 
[numbOfCoord, ~] = maxk(contourPts,50);
bigCidx = []; contourCoordinates = [];
for ii = 3:size(numbOfCoord,2)
    bigCidx = [bigCidx find(M(2,:) == numbOfCoord(ii))];
    contourCoordinates = [contourCoordinates M(:,bigCidx+1:numbOfCoord(ii) + bigCidx)]; 
end
cc = unique(contourCoordinates','rows')';
ccx = cc(1,:); ccy = cc(2,:);

%Check the contour
line(ccx,ccy,'Color','k')
hold off

%find average across contour
patch(ccx,ccy,'red')
%}
%}
