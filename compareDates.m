function laterDate = compareDates(fileName, dateToComp)
%dateToComp should be in format of: '220131'

dateIdx = strfind(fileName, '_22'); dateIdx = dateIdx(1) + 1;
dateIdx2 = strfind(dateToComp, '22'); dateIdx2 = dateIdx2(1) + 1;

if isempty(dateIdx)
    dateIdx = strfind(fileName, '_21'); 
end
if isempty(dateIdx2)
    dateIdx2 = strfind(dateToComp, '21'); 
end

year  = str2double(['20' fileName(dateIdx:dateIdx+1)]); 
month = str2double(fileName(dateIdx+2: dateIdx+ 3));
day   = str2double(fileName(dateIdx+4: dateIdx + 5));

year2  = str2double(['20' fileName(dateIdx2:dateIdx2+1)]); 
month2 = str2double(fileName(dateIdx2+2: dateIdx2+ 3));
day2   = str2double(fileName(dateIdx2+4: dateIdx2+ 5));

fileDate  = datetime(year, month, day); dateTC = datetime(year2,month2,day2);
laterDate = fileDate >= dateTC ;

end