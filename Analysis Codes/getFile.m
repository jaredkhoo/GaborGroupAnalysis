function filepath = getFile(experimentName, subjectName)
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');

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

            % Extract the data from the csv into an array
            filepath = fullfile(folder, string(csvFiles(j,1)));

        end
    end
end

   
        
