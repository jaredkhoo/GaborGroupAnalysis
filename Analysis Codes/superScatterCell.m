function superScatterCell()
    folder = fullfile(pwd, 'Plots', 'SuperScatter');
    mkdir(folder);
    subjectInfo = table2struct(readtable(fullfile(pwd, 'scatterinfo.csv')));
    info_csv = fullfile(pwd, 'protocol_info.csv');
    protocolInfo = table2struct(readtable(info_csv));
    pmodifier = length(protocolInfo);
    h = 1;
    
    while h <= pmodifier
        if isequal(protocolInfo(h).Include, 0)
            protocolInfo(h) = [];
            pmodifier = pmodifier-1;
        else
            h = h+1;
        end
    end
    numgraphs = length(protocolInfo);
    for r=1:numgraphs
        figure(r);
    end
    hold on;
    g = 1;
    smodifier = length(subjectInfo);
    while g <= smodifier
        if isequal(subjectInfo(g).hasData, 0)
            subjectInfo(g) = [];
            smodifier = smodifier-1;
        else
            g = g+1;
        end
    end
    for i=1:length(protocolInfo)
        experimentName = protocolInfo(i).protocolName;
        discreteCol = protocolInfo(i).discreteCol;
        id = protocolInfo(i).id;
        stim = protocolInfo(i).stim;
        abbrev = protocolInfo(i).abbrev;
        for j=1:length(subjectInfo)
            color = str2num(subjectInfo(j).color);
            subjectName = subjectInfo(j).subjectName;
            graphName = subjectInfo(j).graphName;
            thisfile = getFile(experimentName, subjectName);
            raw = table2cell(readtable(thisfile));
            raw(:,1) = [];
            raw = cell2mat(raw); %%% don't need for non-cell succint analysis, aka
            %not character
            data(:,1) = raw(:,2);
            data(:,2) = raw(:,3).*1000;
            data(:,3) = raw(:,1);
            data(isnan(data)) = 0;
            for k = 1:length(data)
                if data(k,3) == 2
                    data(k,1) = data(k,1)*-1;
                end
            end
            data(:,3) = [];
            [nooutliers, theoutliers] = improvedOutliers(data, [], 2.5, 2);
            avgdata = averageData(nooutliers, 1);
            for u=1:length(avgdata)
                measurement = avgdata(u,1);
                idx = find(avgdata(:,1) == measurement);
                avgdata(u,4) = avgdata(idx(1),3);
            end
            figure(i);
            errorbar(avgdata(:,1), avgdata(:,2), avgdata(:,3), 'vertical','.', ...
                'HandleVisibility', 'off', 'Color', color, ...
                'CapSize', 0);
            hold on;
            scatter(avgdata(:,1),avgdata(:,2),30,color,'filled', 'DisplayName', graphName);
            legend('show', 'location', 'northoutside', 'NumColumns', 5);
            for m=1:2  %1 = negative eccentricities, 2 = positive eccentrities
                if m==1
                    b = 1;
                    for n=1:length(avgdata)
                        if avgdata(n,1) <= 0
                            xvals(b,1) = avgdata(n,1);
                            yvals(b,1) = avgdata(n,2);
                            sErrors(b,1) = avgdata(n,4);
                            b = b+1;
                        end
                    end
                    sErrors(sErrors == 0) = 0.000000001;
                    
                    weights = (1./sErrors).^2; %%%%%%
                    f = @(x, xvals, yvals, weights)sum(weights.*((yvals-((xvals.*x(1))+x(2))).^2));
                    fun = @(x)f(x,xvals,yvals,weights);
                    approx = polyfit(xvals, yvals, 1);
                    approxSlope = approx(1); approxInt = approx(2);
                    ms = MultiStart;
                    problem = createOptimProblem('fmincon','x0',approx, ...
                        'objective',fun, ...
                        'lb', [(approxSlope-10), (approxInt-100)],...
                        'ub',[(approxSlope+10), (approxInt+100)]);     %changed here for ms
                    params = run(ms, problem, 50);
                    slope = params(1);
                    intercept = params(2);
                    chisquare = fun(params);
                    reducedchi = chisquare/(length(xvals)-2);

                    figure(i);


                    if xvals(1,1) < 0
                        xfit = linspace(min(data(:,1)), 0);
                    else
                        xfit = linspace(0, max(data(:,1)));
                    end

                    yfit = (xfit*slope)+intercept;

                    space = ' ';
                    negative = '(Negative Eccentricity)';


                    if intercept == 0
                        txt = "%s: y = %4.4fx";

                        plot(xfit, yfit, 'Color', color, 'LineWidth', 1, 'DisplayName', ...
                            sprintf(txt, strcat(experimentName, space, negative), slope));

                    else
                        %txt = "%s: y = %4.4fx + %4.4f";
                        txt = "%s";
%changed here
                        plot(xfit, yfit, 'Color', color, 'LineWidth', .5, 'LineStyle', ':', 'HandleVisibility','off');
                    end

                    xlabel("Retinal Eccentricity (degrees)");
                    ylabel("Reaction Time (ms)"); %changed here for ms
                    title("Reaction time vs. Retinal Eccentricity");


                    xlim([-40 40])
                    ylim([100 900]) %changed here for ms
                    hold on;
                    
                else
                    a = 1;
                    for n=1:length(avgdata)
                        if avgdata(n,1) >= 0
                            xvals(a,1) = avgdata(n,1);
                            yvals(a,1) = avgdata(n,2);
                            sErrors(a,1) = avgdata(n,4);
                            a = a+1;
                        end
                    end
                    sErrors(sErrors == 0) = 0.000000001;
                    weights = (1./sErrors).^2;
                    f = @(x, xvals, yvals, weights)sum(weights.*((yvals-((xvals.*x(1))+x(2))).^2));
                    fun = @(x)f(x,xvals,yvals,weights);
                    approx = polyfit(xvals, yvals, 1);
                    approxSlope = approx(1); approxInt = approx(2);
                    ms = MultiStart;
                    problem = createOptimProblem('fmincon','x0',approx, ...
                        'objective',fun, ...
                        'lb', [(approxSlope-10), (approxInt-100)],...
                        'ub',[(approxSlope+10), (approxInt+100)]);  %changed here for ms
                    params = run(ms, problem, 50);
                    slope = params(1);
                    intercept = params(2);
                    chisquare = fun(params);
                    reducedchi = chisquare/(length(xvals)-2); 

                    figure(i);



                    if xvals(1,1) < 0
                        xfit = linspace(min(data(:,1)), 0);
                    else
                        xfit = linspace(0, max(data(:,1)));
                    end

                    yfit = (xfit*slope)+intercept;

                    positive = '(Positive Eccentricity)';
                    space = ' ';

                    if intercept == 0
                        txt = "%s: y = %4.4fx";

                        plot(xfit, yfit, 'Color', color, 'LineWidth', 1, 'DisplayName', ...
                            sprintf(txt, strcat(experimentName, space, positive), slope));

                    else
                        %txt = "%s: y = %4.4fx + %4.4f";
%%% Changed here
                        plot(xfit, yfit, 'Color', color, 'LineWidth', .5, 'LineStyle', ':', 'HandleVisibility','off');
                    end


                    xlabel("Retinal Eccentricity (degrees)");
                    ylabel("Reaction Time (ms)");   %changed here for ms
                    title("Reaction time vs. Retinal Eccentricity");


                    xlim([-40 40])
                    ylim([100 900])
                    hold on;
                end
            end
        clear data
        end
        finalfig = strcat(abbrev,'_SuperScatter_', stim, '.png');
        saveas(figure(i), fullfile(folder, finalfig));
            close(figure(i));
    end

end



