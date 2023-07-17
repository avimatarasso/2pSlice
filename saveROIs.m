%% Save ROIs and our images from each session.

clear all

addpath('G:\code\2pSliceAnalysis\oir2stdData-master');
%addpath('C:\Users\avima\OneDrive\Documents\GitHub\2pSlice');
addpath('G:\code\2pSliceAnalysis\');
pathToFile = 'G:\2P-slice exps\Sensor Control Exps\GRABDA\216740-1\stim';%'F:\2P-slice exps\NE DA proj\dlight\washes';%'F:\2P-slice exps\022222\';%'F:\2P-slice exps\022222\105304-3_S1_dL_CA1_20hz_30s_15sdelay.oir';
pathToFile = 'G:\2P-slice exps\Sensor Control Exps\GRABDA\216740-1';
cd(pathToFile)
filePattern = ['216740'];
workFiles = dir([pathToFile '\*' filePattern '*.oir']);
stimON = 0;
skip = 0;
multiple = 1; 
justSave = 0; %if you dont want to compare
ROIdir = 'saved ROIs';
if ~exist(ROIdir,'dir')
    mkdir(ROIdir)
end

oirFiles = dir('*.oir');
% Goal: create a mask for each slice, then check how the mask looks within
% one session and confirm it looks good, and save that mask under each 

ptr = 1;
for fileN=1:length(oirFiles)
    currentOir = oirFiles(fileN).name; oirSlice = strsplit(currentOir,'_');
    roiFilename = [currentOir(1:end-4) '_roi.mat'];
    
    % uncomment below to skip already made rois
    if exist([pwd '/saved ROIs/' roiFilename ],'file') && skip
        continue
    end
    
    %find where slice naming ends, may need to CUSTOMIZE
    delimIdx = strfind(currentOir,'_'); slcStrEnd = delimIdx(2)-1; 
    currentSlice = currentOir(1:slcStrEnd);
        
    %[~,stdData]=oir2stdData(currentOir);
    %a = squeeze(stdData(1).Image{1});
    if stimON
        [sessionImgs, stdData, stimTime, delay, freq] = initialize2p(currentOir,stimON); 
        stimTime = details.stimTime;
        delay = details.delay;
        freq = details.freq;
        cd('saved ROIs') 
        save(roiFilename,'sessionImgs', 'stdData','stimTime', 'delay', 'freq','-v7.3')
        cd('..')
    else
         [sessionImgs, stdData] = initialize2p(currentOir,stimON);
         cd('saved ROIs')
         save(roiFilename,'sessionImgs', 'stdData','-v7.3')
        cd('..')
    end
    
    disp(['Now creating an ROI for: ' currentSlice])
    disp(['With file: ' currentOir])


    roiFilenameGen = dir(['saved ROIs\' currentSlice '*.mat']);
    load([pwd '/saved ROIs/' roiFilename])
    if ~exist('ROImask','var')
        % if there is no roi files with this slice in 'saved ROIs', create a
        % mask, if the roi file exists (else) then check the roi file matches
        meanIMG = mean(sessionImgs(:,:,[1 end]),3);
        cd('saved ROIs')
        ROImask = createROImask(meanIMG); 
        ROImask = compareROI(meanIMG,ROImask,currentOir,justSave);        
        save(roiFilename, 'ROImask','meanIMG','-v7.3','-append')
    else
        % ask does the ROI match the image?
        % if ROI does not match, then create a new mask, and if it does
        % match, then use the same mask and save a new roi file for that slice!
        cd('saved ROIs')
        if multiple
            tmpFilename = [currentOir(1:end-4) '_roi.mat']; %'_000' num2str(ptr) ?
        else
            tmpFilename = [currentOir(1:end-4) '_roi.mat'];
        end
        load(tmpFilename)
        %overwrite meanIMG of the original roiFileName for currentOir
        meanIMG = mean(sessionImgs(:,:,[1 end]),3);
        ROImask = compareROI(meanIMG,ROImask,currentOir,justSave);       
        %if multiple && fileN == 1
        %    ROImaskFull = ROImask;
        %elseif multiple && fileN ~=1
        %    ROImaskFull = ROImaskFull + ROImask;
        %end
        save(tmpFilename,'ROImask','meanIMG','-v7.3','-append')
    end
    cd('..')    
end
%{
quest    = 'Do you want another ROI?';
dlgtitle = 'Multiple ROIs?';
btn1     = 'Yes, more ROIs!'; 
btn2     = 'No, we captured every ROI'; defbtn = btn2; 
answer  = questdlg(quest,dlgtitle,btn1,btn2,defbtn); 

%make multiple == 0 when done!
%confirm with a comparison image of all ROIs
if strcmp(answer, btn2)
    switchBtn = 0;
    
    figure(300)
    fig.Position = [100 100 1000 800];
    % Put both images into subplot
    subplot(1,2,1)
        imshow(meanImg, [], 'Colormap', gray(256));
    subplot(1,2,2)
    ha = imshow(ROImaskFull, [], 'Colormap', gray(256));
    hold on;
    hb = imshow(ROImaskFull, [], 'Colormap', gray(256));
    sgtitle(strrep(titleFile,'_',' '))
    % Set opacity/transparency to something less than 1 (alpha).  
    % 1 is the default and it means the last image is opaque and the image below can't be seen.
    hb.AlphaData = 0.1;
    
    quest    = 'Did you capture all ROIs?';
    dlgtitle = 'Multiple ROIs DONE?!';
    btn1     = 'Yes, we"re done!'; 
    btn2     = 'No, we need to do more'; defbtn = btn2; 
    answer  = questdlg(quest,dlgtitle,btn1,btn2,defbtn); 
    if strcmp(answer,btn2)
        switchBtn = 1;
    end
end
%}
%ptr = ptr+1;
%end