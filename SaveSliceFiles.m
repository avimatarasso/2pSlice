
%% Save 2p files
addpath('F:\code\2pSliceAnalysis\oir2stdData-master');
addpath('C:\Users\avima\OneDrive\Documents\GitHub\2pSlice');
pathToFile = 'F:\2P-slice exps\NE DA proj\GRABNE\washes';%'F:\2P-slice exps\022222\';%'F:\2P-slice exps\022222\105304-3_S1_dL_CA1_20hz_30s_15sdelay.oir';
cd(pathToFile)
filePattern = 'CA1';
workFiles = dir([pathToFile '\*' filePattern '*.oir']);
stimON = 0;

for fileN = 1:length(workFiles)
    fileName = workFiles(fileN).name; 
    fileNameMat = [fileName(1:end-4) '.mat'];
    if ~exist(fileNameMat)
    % extract the metadata and the files
        if stimON
            [fileName, CRGB, stimTime, delay] = initialize2p(workFiles, fileN,stimON);
            save(fileNameMat, 'CRGB', 'stimTime', 'delay')
        else
            [fileName, CRGB] = initialize2p(workFiles, fileN,stimON);
            save(fileNameMat, 'CRGB', '-v7.3')            
        end
    end
end
