function createAllCellType(protocolInfo)
    addpath(fullfile(pwd));
    sub = getSubjects();
    filestruct = struct('name', NaN, 'filepath', NaN, 'exp', NaN);
    for c=1:length(protocolInfo)
        counter = 1;
        experimentName = protocolInfo(c).protocolName;
        for a=2:length(sub)
            warning('off','MATLAB:table:ModifiedAndSavedVarnames');
            subjectName = sub(a).name;

            % Creates filepaths to all the subfolders within each base folder
            folderPaths = [];

            baseFolder = string(fullfile(pwd, 'Data'));

            folders = dir(baseFolder);
            folderNames={folders(:).name}';
            for j=1:size(folderNames,1)
                folderName = string(folderNames{j,1});

                if(startsWith(folderName,"."))
                    continue; % Skips folders such as .DS_STORE
                end

                % If a subject was specified, skip any folders which do not
                % match the subject code
                if(~strcmp(subjectName,'All'))
                    if(~strcmp(subjectName,folderName))
                        continue;
                    end
                end

                % Create the new file path and append it to the array
                folderPath = string(fullfile(baseFolder, folderName));
                folderPaths = [folderPaths folderPath];  
            end

            rawData = []; data = []; outliers = [];
            for i=1:length(folderPaths)
                folder = folderPaths(i);

                % Get the names of all the csv files within the folder
                files = dir(folder);
                fileNames={files(:).name}';
                csvFiles=fileNames(endsWith(fileNames,'.csv'));

                % Loop over every csv file, extracting data from files which match
                % the search conditions
                for j=1:size(csvFiles,1)
                    file = char(csvFiles{j,1});

                    % csv file names follow this convention:
                    % SUBJECT_DATE_PROTOCOLNAME.csv
                    underscores = find(file == '_'); % Find indices of underscores
                    period = find(file == '.'); % Find the period

                    % Extract the portion of the csv file name between the last
                    % underscore and the period
                    protocolName = string(extractBetween(file, (underscores(end)+1), (period-1)));
                    if(~strcmp(protocolName, experimentName))
                        continue; % Skip csv files from other protocols
                    end

                    thisData = []; thisFitData = []; theseOutliers = [];

                    % Extract the data from the csv into an array
                    
                    filename = fullfile(folder, string(csvFiles(j,1)));
                    filestruct(counter).name = sub(a).name;
                    filestruct(counter).filepath = filename;
                    filestruct(counter).exp = experimentName;
                    counter = counter+1;
                end
            end
        end
        folder = fullfile(pwd, 'Data', 'All');
        mkdir(folder);
        file = fullfile(pwd, 'Data', 'All', strcat('All_', experimentName, '.csv'));
        for k=1:length(filestruct)
            data = table2cell(readtable(filestruct(k).filepath)); %changed from table2array
            if(exist(file, 'file') ~= 2)
                writecell(data, file); %
            else
                writecell(data, file, 'WriteMode', 'Append');
            end
        end
        w = 1;
        x = length(filestruct);
        while w <= x
            filestruct(x) = [];
            x = x-1;
        end
        
    end
end
        
    