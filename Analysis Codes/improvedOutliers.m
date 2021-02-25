function [finaldata, outliers] = improvedOutliers(data, outliers, cutOff, column)

    sorteddata = sortrows(data, 1);
    finaldata = [];
    tempdata = [];
    a = unique(sorteddata(:,1));
    lenfinal = length(finaldata)+1;
    for j = 1:length(a)
        k = 1;
        count = 1;
        while k <= length(sorteddata)
            if sorteddata(k,1) == a(j,1)
                tempdata(count,:) = sorteddata(k,:);
                count = count+1;
            end
            k = k + 1;
        end
        z = zscore(tempdata(:,column));
        tempdata(:,column+1) = z(:,1);
        l = 1;
        while(l <= length(tempdata))
            if (abs(tempdata(l,column+1))) > cutOff
                outliers = [outliers; tempdata(l,:)];
                tempdata(l,:) = [];
            end
            l = l+1;
        end
        tempdata(:,column+1) = [];
        for m = 1:length(tempdata)
            finaldata(lenfinal,:) = tempdata(m,:);
            lenfinal = lenfinal + 1;
        end
        tempdata = [];
    end
    if length(sorteddata) ~= length(finaldata)
        [finaldata, outliers] = improvedOutliers(finaldata, outliers, cutOff, column);
    end
end

            
        