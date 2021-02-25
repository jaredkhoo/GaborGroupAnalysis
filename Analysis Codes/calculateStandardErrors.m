% calculateStandardErrors
%
% Estimates the standard errors of each observation for the Chi^2
% minimization as the standard error of the distribution of all
% observations ever recorded at that discrete measurement, across
% all subjects. Takes a struct containing information about the protocol,
% a boolean indicating whether small eccentricity observations should be
% excluded from crowded center, and the data matrix.
function outputData = calculateStandardErrors(experimentName, discreteCol, outputData)
    
    % Extract all of the data ever recorded for the current protocol from
    % the csv files
    [~,data,~] = readCsv(experimentName,'All','both');
    
    % Convert the data back into the linear scale (from dimensionless y/x),
    % and average all observations made at each discrete measurement point
    data(:,2) = data(:,2).*data(:,1);
    avgData = averageData(data,discreteCol);
    
    % For each measurement in the data matrix, find it's matching entry in
    % the avgData matrix and copy the standard error from that matrix
    for i=1:length(outputData)
        measurement = outputData(i,discreteCol);
        idx = find(avgData(:,discreteCol) == measurement);
        outputData(i,4) = avgData(idx(1),3);
    end
    
end