function [CRGB, stdData, details] = initialize2p(fileName,stimON)
strfind(fileName,'s_');
splitStr = regexp(fileName,'_','split');
len  = length(splitStr);

if stimON
freq = splitStr(contains(splitStr,'hz')); details.freq= freq{1};
delayStr = splitStr(contains(splitStr,'sdelay')); delayStr = delayStr{1};
stimTimeStr = splitStr(find(contains(splitStr,'sdelay'))-1); stimTimeStr = stimTimeStr{1};
stimTime = stimTimeStr(1:strfind(stimTimeStr,'s')-1);
details.stimTime = str2double(stimTime);
delay = delayStr(1:strfind(delayStr,'sdelay')-1); %make it in seconds!
details.delay = str2double(delay);    
end

[~,stdData]=oir2stdData(fileName);
a = squeeze(stdData(1).Image{1});
% convert uint8 to RGB to visualize
CRGB = zeros(size(a,1),size(a,2),size(a,3));
for imgNumb = 1:size(a,3)
    C = mat2gray(a(:,:,imgNumb),[0 4096]);
    CRGB(:,:,imgNumb) = C;
end


end