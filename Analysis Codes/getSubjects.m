function subjects = getSubjects()

    subjects = struct('name', 'All', 'includeAll', true);

    folders = dir(fullfile(pwd, 'Data'));

    folderNames={folders(:).name}';

    for j=1:size(folderNames,1)
        folderName = string(folderNames{j,1});

        if(startsWith(folderName,"."))
            continue; % Skips folders such as .DS_STORE
        end
        subjects(end+1) = struct('name', folderName, 'includeAll', true);
    end

end