% avgData
%
% Takes a vertical matrix of data and the discrete column (usually the
% independent variable/measurement) as arguments. Calculates the mean and
% the standard error of the distribution of observations at each unique
% measurement. Returns a shortened matrix, with the standard errors in the
% third column.
function avgData = averageData(data, discreteCol)
    
    % If discreteCol = 2, avgCol = 1. If discreteCol = 1, avgCol = 2
    avgCol = 3-discreteCol;
    
    % Round values in the discrete column to two decimal places to correct
    % any change in value induced by the (y./x).*x transformation
    %if discreteCol == 2
        %data(:,2) = round(data(:,2),2);
    %end
    
    % Get all the unique values in the discrete column, loop through
    avgData(:,discreteCol) = unique(data(:,discreteCol));
    for i = 1:size(avgData,1)
        
        % Find all indices in the raw data matrix where the current value
        % can be found
        indices = find(data(:,discreteCol) == avgData(i,discreteCol));
        values = zeros(1,length(indices));
        
        % Extract all of these values from the raw matrix into an array
        for j = 1:length(indices)
            index = indices(j,1);
            values(j) = data(index,avgCol);
        end
        
        % Average the values and calculate the standard error
        avgData(i,avgCol) = mean(values);
        avgData(i,3) = std(values)/sqrt(length(values));
    end
end