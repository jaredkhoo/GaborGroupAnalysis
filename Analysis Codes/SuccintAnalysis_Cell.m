addpath(fullfile(pwd));

subjects = getSubjects();

csv = fullfile(pwd, 'protocol_info.csv');
pInfo = table2struct(readtable(csv));
pmodifier = length(pInfo);
h = 1;

while h <= pmodifier
    if isequal(pInfo(h).Include, 0)
        pInfo(h) = [];
        pmodifier = pmodifier-1;
    else
        h = h+1;
    end
end
stimname = pInfo(1).stim;
paramfolder = fullfile(pwd, 'Parameters');
mkdir(paramfolder);
fileName = fullfile(pwd, 'Parameters', strcat(stimname, '_chi_statistics.csv'));


for i=1:length(subjects)
    subject = subjects(i);
    analysisoutput = analyzeIndCell(subject.name);
    output = struct2table(analysisoutput);
    if i == 1    
        if(exist(fileName, 'file') ~= 2)
            writetable(output, fileName, 'WriteRowNames', true);
        else
            writetable(output, fileName, 'WriteRowNames', true, ...
                'WriteMode', 'overwrite');
        end
    else
        writetable(output, fileName, 'WriteRowNames', false, ...
            'WriteMode', 'Append');      
    end
end
superScatterCell();