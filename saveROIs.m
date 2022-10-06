%% Save ROIs and our images from each session.

addpath('/Users/catazamorano/Desktop/code/2pSlice-main');
%addpath('C:\Users\avima\OneDrive\Documents\GitHub\2pSlice');
addpath('/Volumes/Cat Seagate/mlight/2 photon');
pathToFile = '/Volumes/Cat Seagate/mlight/2 photon/151649-6';
cd(pathToFile)
filePattern = 'DAMGO';
workFiles = dir([pathToFile '\*' filePattern '*.oir']);
stimON = 0;

% If you want to skip already made slices make skip == 1
skip = 1;

switchBtn = 1; multiple = 0; 
justSave = 0; %if you dont want to compare
ROIdir = 'ROIdir';
if ~exist(ROIdir,'dir')
    mkdir ROIdir
end

oirFiles = dir('*.oir');
% Goal: create a mask for each slice, then check how the mask looks within
% one session and confirm it looks good, and save that mask under each 

%while switchBtn == 1
    
for fileN=1:length(oirFiles)
    currentOir = oirFiles(fileN).name; oirSlice = strsplit(currentOir,'_');
    roiFilename = [currentOir(1:end-4) '_roi.mat'];
    
    % uncomment below to skip already made rois
    if exist([pwd '/ROIdir/' roiFilename ],'file') && skip
        continue
    end
    
    %find where slice naming ends, may need to CUSTOMIZE
    delimIdx = strfind(currentOir,'_'); slcStrEnd = delimIdx(2)-1; 
    currentSlice = currentOir(1:slcStrEnd);
        
    [~,stdData]=oir2stdData(currentOir);
    a = squeeze(stdData(1).Image{1});
    if stimON
        [sessionImgs, stdData, stimTime, delay, freq] = initialize2p(currentOir,stimON);        
        cd('ROIdir') 
        save(roiFilename,'sessionImgs', 'stdData','stimTime', 'delay', 'freq','-v7.3')
        cd('..')
    else
         [sessionImgs, stdData,~,~,~] = initialize2p(currentOir,stimON);
         cd('ROIdir')
         save(roiFilename,'sessionImgs', 'stdData','-v7.3')
        cd('..')
    end
    
    disp(['Now creating an ROI for: ' currentSlice])
    disp(['With file: ' currentOir])


    roiFilenameGen = dir(['ROIdir\' currentSlice '*.mat']);
    load([pwd '/ROIdir/' roiFilename])
    if ~exist('ROImask','var')
        % if there is no roi files with this slice in 'ROIdir', create a
        % mask, if the roi file exists (else) then check the roi file matches
        meanIMG = mean(sessionImgs(:,:,[1 end]),3);
        cd('ROIdir')
        ROImask = createROImask(meanIMG); 
        ROImask = compareROI(meanIMG,ROImask,currentOir,justSave);        
        save(roiFilename, 'ROImask','meanIMG','-v7.3','-append')
    else
        % ask does the ROI match the image?
        % if ROI does not match, then create a new mask, and if it does
        % match, then use the same mask and save a new roi file for that slice!
        cd('ROIdir')
        if multiple
            tmpFilename = [currentOir(1:end-4) '_000' num2str(ptr) '_roi.mat'];
        else
            tmpFilename = [currentOir(1:end-4) '_roi.mat'];
        end
        load(tmpFilename)
        %overwrite meanIMG of the original roiFileName for currentOir
        meanIMG = mean(sessionImgs(:,:,[1 end]),3);
        ROImask = compareROI(meanIMG,ROImask,currentOir,justSave);       
        if multiple && fileN == 1
            ROImaskFull = ROImask;
        elseif multiple && fileN ~=1
            ROImaskFull = ROImaskFull + ROImask;
        end
        save(tmpFilename,'ROImask','meanIMG','-v7.3','-append')
    end
    cd ..    
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