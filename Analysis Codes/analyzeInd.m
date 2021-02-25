function output = analyzeInd(subjectName)
    folder = fullfile(pwd, 'Plots', subjectName);
    mkdir(folder);
    
    info_csv = fullfile(pwd, 'protocol_info.csv');
    protocolInfo = table2struct(readtable(info_csv));
    
    template = fullfile(pwd, 'param_template.csv');
    output = table2struct(readtable(template));
    output.name = subjectName;

    hold on;
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

    numexp = length(protocolInfo);
    numgraphs = (numexp.*2)+1;
    for g=1:numgraphs
        figure(g);
    end
    hold on;
    if strcmp(subjectName, 'All')
        createAll(protocolInfo);
    end
    for i=1:length(protocolInfo)
        color = str2num(protocolInfo(i).color);
        experimentName = protocolInfo(i).protocolName;
        discreteCol = protocolInfo(i).discreteCol;
        id = protocolInfo(i).id;
        stim = protocolInfo(i).stim;
        abbrev = protocolInfo(i).abbrev;
        if strcmp(subjectName, 'All')
            thisfile = getAllFile(experimentName, subjectName);
        else
            thisfile = getFile(experimentName, subjectName);
        end
        raw = table2array(readtable(thisfile));
        data(:,1) = raw(:,3);
        data(:,2) = raw(:,4).*1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        data(:,3) = raw(:,2);
        data(isnan(data)) = 0;
        for j = 1:length(data)
            if data(j,3) == 2
                data(j,1) = data(j,1)*-1;
            end
        end
        data(:,3) = [];
        
        [nooutliers, theoutliers] = improvedOutliers(data, [], 2.5, 2);
        avgdata = averageData(nooutliers, 1);
        
        for k=1:length(avgdata)
            measurement = avgdata(k,1);
            idx = find(avgdata(:,1) == measurement);
            avgdata(k,4) = avgdata(idx(1),3);
        end
        
        figure(1);
        errorbar(avgdata(:,1), avgdata(:,2), avgdata(:,3), 'vertical','.', ...
            'HandleVisibility', 'off', 'Color', color, ...
            'CapSize', 0);

        hold on;
        scatter(avgdata(:,1),avgdata(:,2),30,color,'filled', 'DisplayName', abbrev);
        legend('show', 'location', 'eastoutside');
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
                
                weights = (1./sErrors).^2;
                %f = @(x, xPoints, yPoints, w)sum(w.*((yPoints-((xPoints.*x(1))+x(2))).^2));
                %optFun = @(x)f(x, xvals, yvals, weights);
                %ms = MultiStart;
                %OLSFit = polyfit(xvals, yvals, 1);
                corCoef = corrcoef(xvals,  yvals);
                R = corCoef(1,2);
                %guessParams = [OLSFit(1), OLSFit(2)];
                %problem = createOptimProblem('fmincon', 'x0', guessParams, ...
                %    'objective', optFun, 'lb' , [-15, 100], 'ub', [15, 700] );
                %params = run(ms,problem,25);
                %slope = params(1);
                %intercept = params(2);
                %chiSquareFit.slope = slope;
                %chiSquareFit.intercept = intercept;
                %chi2Val = optFun(params);
                %chisquare = optFun(params);
                
                %weights = 1./sErrors;
                %weights = weights.^2;
                f = @(x, xvals, yvals, weights)sum(weights.*((yvals-((xvals.*x(1))+x(2))).^2));
                fun = @(x)f(x,xvals,yvals,weights);
                approx = polyfit(xvals, yvals, 1);
                approxSlope = approx(1); approxInt = approx(2);
                ms = MultiStart;
                problem = createOptimProblem('fmincon','x0',approx, ...
                    'objective',fun, ...
                    'lb', [(approxSlope-10), (approxInt-100)],...
                    'ub',[(approxSlope+10), (approxInt+100)]);      %changed here for ms
                params = run(ms, problem, 50);
                slope = params(1);
                intercept = params(2);
                chisquare = fun(params);
                reducedchi = chisquare/(length(xvals)-2);
                
                figure(1);
                
               

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

                    plot(xfit, yfit, 'Color', color, 'LineWidth', 1, 'HandleVisibility','off');
                end

                xlabel("Retinal Eccentricity (degrees)");
                ylabel("Reaction Time (ms)"); %changed here for ms
                title("Reaction time vs. Retinal Eccentricity");


                xlim([-40 40])
                ylim([100 900]) %changed here for ms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                hold on;
                
                target = chisquare+1;
                f = @(x, intercept, xvals, yvals, weights)abs(target-sum(weights.*((yvals-((xvals.*x)+intercept)).^2)));
                fun = @(x)f(x, intercept, xvals, yvals, weights);

                if slope < 0
                    negSlopeError = fminbnd(fun, (slope*1.5), slope);   %might have to change back to 0.5
                    posSlopeError = fminbnd(fun, slope, (slope*0.5));   %and 1.5 respectively when doing positive
                else
                    negSlopeError = fminbnd(fun, (slope*0.5), slope);   %might have to change back to 0.5
                    posSlopeError = fminbnd(fun, slope, (slope*1.5));%values
                end
                f = @(slope, x, xvals, yvals, weights)abs(target-sum(weights.*((yvals-((xvals.*slope)+x)).^2)));
                fun = @(x)f(slope, x, xvals, yvals, weights);

                negIntError = fminbnd(fun, -1, intercept);
                posIntError = fminbnd(fun, intercept, intercept + 1);  %changed here for ms

                f = @(slope, intercept, xvals, yvals, weights)sum(weights.*((yvals-((xvals.*slope)+intercept)).^2));
                fun = @(slope, intercept)f(slope, intercept, xvals, yvals, weights);

    %check code if need to add back

                slope_range = [negSlopeError posSlopeError];
                int_range = [negIntError posIntError];

                complete = false;
                while(~complete)
                    if slope_range(2) > slope_range(1)
                        slope_evals = linspace(slope_range(1), slope_range(2)); 
                    else
                        slope_evals = linspace(slope_range(2), slope_range(1));%change this back to 1, 2 for positive
                    end
                    int_evals = linspace(int_range(1), int_range(2)); %changed here for ms (the slope range)

                    [slope_grid, int_grid] = meshgrid(slope_evals, int_evals);

                    chi_grid = cellfun(fun, num2cell(slope_grid), num2cell(int_grid));

                    [slope_range,int_range,complete] = checkGrid(chi_grid,slope_range, ...
                        int_range,(chisquare+4));
                end
                negfigure = i.*2;
                figure(negfigure);
                ax = axes('Parent', figure(negfigure));

                h = surf(slope_grid, int_grid, chi_grid, 'Parent', ax,...
                    'edgecolor', 'none');

                hold on;

                contour3(slope_grid, int_grid, chi_grid, ...
                    [(chisquare +1) (chisquare +1)], ...
                    'ShowText', 'off', 'LineColor', 'w', 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');

                view(ax, [0, 90]);
                colormap(jet);

                chisqtext = "\chi^{2}"; 

                ccb = colorbar;
                ccb.Label.String = chisqtext;
                ccb.Label.FontSize = 12;
                ccb.Label.Rotation = 0;
                pos = get(ccb.Label, 'Position');
                ccb.Label.Position = [pos(1)*1.3 pos(2) 0];

                xlabel("Slope");
                ylabel("Intercept");
                txt = "%s vs. Slope and Intercept Parameters";
                title(sprintf(txt, chisqtext));
                if slope_range(1)~=slope_range(2)
                    xlim(slope_range);
                end
                ylim(int_range);
        
                output.(strcat('slope_', id, '_neg')) = slope;
                output.(strcat('intercept_', id, '_neg')) = intercept;
                output.(strcat('chi_', id, '_neg')) = chisquare;
                output.(strcat('reducedchi_', id, '_neg')) = reducedchi;
                output.(strcat('r_', id, '_neg')) = R;
                output.(strcat('negslopeerror_', id, '_neg')) = negSlopeError;
                output.(strcat('posslopeerror_', id, '_neg')) = posSlopeError;
                output.(strcat('neginterror_', id, '_neg')) = negIntError;
                output.(strcat('posinterror_', id, '_neg')) = posIntError;
                saveas(figure(negfigure), fullfile(folder, ... 
                    strcat(subjectName, '_', experimentName, '.png')));
                close(figure(negfigure));

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
                %f = @(x, xPoints, yPoints, w)sum(w.*((yPoints-((xPoints.*x(1))+x(2))).^2));
                %optFun = @(x)f(x, xvals, yvals, weights);
                %ms = MultiStart;
                %OLSFit = polyfit(xvals, yvals, 1);
                corCoef = corrcoef(xvals,  yvals);
                R = corCoef(1,2);
                %guessParams = [OLSFit(1), OLSFit(2)];
                %problem = createOptimProblem('fmincon', 'x0', guessParams, ...
                 %   'objective', optFun, 'lb' , [-15, 100], 'ub', [15, 700] );
                %params = run(ms,problem,25);
                %slope = params(1);
                %intercept = params(2);
                %chiSquareFit.slope = slope;
                %chiSquareFit.intercept = intercept;
                %chi2Val = optFun(params);
                %chisquare = optFun(params);
                
                %weights = 1./sErrors;
                %weights = weights.^2;
                f = @(x, xvals, yvals, weights)sum(weights.*((yvals-((xvals.*x(1))+x(2))).^2));
                fun = @(x)f(x,xvals,yvals,weights);
                approx = polyfit(xvals, yvals, 1);
                approxSlope = approx(1); approxInt = approx(2);
                ms = MultiStart;
                problem = createOptimProblem('fmincon','x0',approx, ...
                    'objective',fun, ...
                    'lb', [(approxSlope-10), (approxInt-100)],...
                    'ub',[(approxSlope+10), (approxInt+100)]);      %changed here for ms
                params = run(ms, problem, 50);
                slope = params(1);
                intercept = params(2);
                chisquare = fun(params);
                reducedchi = chisquare/(length(xvals)-2);
                
                figure(1);



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

                    plot(xfit, yfit, 'Color', color, 'LineWidth', 1,'HandleVisibility','off');
                end


                xlabel("Retinal Eccentricity (degrees)");
                ylabel("Reaction Time (ms)");   %changed here for ms
                title("Reaction time vs. Retinal Eccentricity");


                xlim([-40 40])
                ylim([100 900]) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hold on;
                
                target = chisquare+1;
                f = @(x, intercept, xvals, yvals, weights)abs(target-sum(weights.*((yvals-((xvals.*x)+intercept)).^2)));
                fun = @(x)f(x, intercept, xvals, yvals, weights);

                if slope < 0
                    negSlopeError = fminbnd(fun, (slope*1.5), slope);   %might have to change back to 0.5
                    posSlopeError = fminbnd(fun, slope, (slope*0.5));   %and 1.5 respectively when doing positive
                else
                    negSlopeError = fminbnd(fun, (slope*0.5), slope);   %might have to change back to 0.5
                    posSlopeError = fminbnd(fun, slope, (slope*1.5));%values
                end

                f = @(slope, x, xvals, yvals, weights)abs(target-sum(weights.*((yvals-((xvals.*slope)+x)).^2)));
                fun = @(x)f(slope, x, xvals, yvals, weights);

                negIntError = fminbnd(fun, -1, intercept);
                posIntError = fminbnd(fun, intercept, intercept + 1);

                f = @(slope, intercept, xvals, yvals, weights)sum(weights.*((yvals-((xvals.*slope)+intercept)).^2));
                fun = @(slope, intercept)f(slope, intercept, xvals, yvals, weights);

    %check code if need to add back

                slope_range = [negSlopeError posSlopeError];
                int_range = [negIntError posIntError];

                complete = false;
                while(~complete)
                    if slope_range(2) > slope_range(1)
                        slope_evals = linspace(slope_range(1), slope_range(2));  %change this back to 1, 2 for positive
                    else
                        slope_evals = linspace(slope_range(1), slope_range(2));
                    end
                    int_evals = linspace(int_range(1), int_range(2));

                    [slope_grid, int_grid] = meshgrid(slope_evals, int_evals);

                    chi_grid = cellfun(fun, num2cell(slope_grid), num2cell(int_grid));

                    [slope_range,int_range,complete] = checkGrid(chi_grid,slope_range, ...
                        int_range,(chisquare+4));
                end
                posfigure = ((i.*2)+1);
                figure(posfigure);
                ax = axes('Parent', figure(posfigure));
                
                h = surf(slope_grid, int_grid, chi_grid, 'Parent', ax,...
                    'edgecolor', 'none');

                hold on;

                contour3(slope_grid, int_grid, chi_grid, ...
                    [(chisquare +1) (chisquare +1)], ...
                    'ShowText', 'off', 'LineColor', 'w', 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');

                view(ax, [0, 90]);
                colormap(jet);

                chisqtext = "\chi^{2}"; 

                ccb = colorbar;
                ccb.Label.String = chisqtext;
                ccb.Label.FontSize = 12;
                ccb.Label.Rotation = 0;
                pos = get(ccb.Label, 'Position');
                ccb.Label.Position = [pos(1)*1.3 pos(2) 0];

                xlabel("Slope");
                ylabel("Intercept");
                txt = "%s vs. Slope and Intercept Parameters";
                title(sprintf(txt, chisqtext));
                if slope_range(1)~=slope_range(2)
                    xlim(slope_range);
                end
                ylim(int_range);
                
                output.(strcat('slope_', id, '_pos')) = slope;
                output.(strcat('intercept_', id, '_pos')) = intercept;
                output.(strcat('chi_', id, '_pos')) = chisquare;
                output.(strcat('reducedchi_', id, '_pos')) = reducedchi;
                output.(strcat('r_', id, '_pos')) = R;
                output.(strcat('negslopeerror_', id, '_pos')) = negSlopeError;
                output.(strcat('posslopeerror_', id, '_pos')) = posSlopeError;
                output.(strcat('neginterror_', id, '_pos')) = negIntError;
                output.(strcat('posinterror_', id, '_pos')) = posIntError;
                saveas(figure(posfigure), fullfile(folder, ... 
                    strcat(subjectName, '_', experimentName, '.png')));
                close(figure(posfigure));

            end

        end
    clear data
    end
    finalfig = strcat(subjectName,'_chisquarelinearfit_', stim, '.png');
    saveas(figure(1), fullfile(folder, finalfig));
        close(figure(1));

end
        

        
        
        
        
        
    