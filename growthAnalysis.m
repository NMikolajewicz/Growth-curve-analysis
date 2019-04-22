%% Growth Curve Analysis
% by N Mikolajewicz (11.02.19)

clear all; close all;

%% specify inputs and analysis parameters
file = 'growthData.xlsx';           % name of excel file where input data is stored
sheet = 'example';                  % name of excel sheet

interpInterval = 5;                 % specify interval between  days
%(e.g., for days 0,1,2,3, the interval is 1. For 0,2,4,6, the interval is 2.)

% specify range of dates to analyze (e.g., from 0 to 140 days)
startDate = 0;                      % e.g., day 0
endDate = 140;                      % e.g., day 140

weightScale = 2;
% 1 = day 0 is baseline weight
% 2 = substract baseline weight from all weight
% 3 = substract min weight from all weights
% 4 = percent different from baseline
% 5 = percent different from minimum weight

startTime = 2;
% 1: earliest date at which measurement is recorded is used as t = 0;
% 2: time at which smallest weight is recorded is used as t = 0;

saveResults = false; % true or false; save results to excel sheet


%% do not modify

% define range of data that is required to be included in analysis
xx = [startDate:interpInterval:endDate]; %

%import data from excel spreadsheet. Refer to template for spreadsheet format.
data = importData(file, sheet); 

% recode categorical variables to numerical variables
categoricalVar = {'condition'}; 
[data, codingLegend] = cat2num(data, categoricalVar);
p=numSubplots(length(unique([data.condition])));

% convert dates to serial date numbers
dateVar = {'startDate', 'weightDate', 'birthDate'}; 
[data] = date2time(data, dateVar);

% idetnfiy unique mouse IDs
uMice = unique([data.ID]);

figure(10);
for i = 1:length(uMice)
    
    % organize data in data structure
    pData(i).ID = uMice(i);
    pData(i).condition = unique([data([data.ID] == uMice(i)).condition]);
    if length([ pData(i).condition]) > 1;  pData(i).condition = nan(); end
    pData(i).startDate = unique([data([data.ID] == uMice(i)).startDate]);
    
    % if no start date is available, use first date of weighing as start date.
    if isnan(pData(i).startDate)
        pData(i).startDate = min([data([data.ID] == uMice(i)).weightDate]);
        if isnan(pData(i).startDate)
            msgbox('check your inputs, there seems to be a problem with the start date variable.');
        end
    end
    
    pData(i).birthDate = unique([data([data.ID] == uMice(i)).birthDate]);
    pData(i).weightDate = [data([data.ID] == uMice(i)).weightDate];
    pData(i).weight = [data([data.ID] == uMice(i)).weight];
    
    % clean up data
    x = [pData(i).weightDate];
    y = [pData(i).weight];
    [x,y] = prepareCurveData(x,y);
  
    % assign start time
    if startTime == 1 % 1: earliest date at which measurement is recorded is used as x = 0;
        x = x-min(x);
    elseif startTime ==2 % 2: time at which smallest weight is recorded is used as x = 0;
        [ymin, yminInd] = min(y); % find min weight, and use index to find time of min weight
        x = x-x(yminInd);
    end
    
    y = y(x>=0); % ensure no date for negative dates (x<0) is included
    x = x(x>=0);
    baselineWeight(i) = y(1);

    % transform weights (as specified)
    if weightScale == 1
    elseif weightScale == 2
        y = y-baselineWeight(i); % subtract baseline weight from all weights
    elseif weightScale == 3
        y = y-min(y); %  % subtract minimum weight from all weights
    elseif weightScale == 4
        y = 100*(y-baselineWeight(i))/baselineWeight(i);
    elseif weightScale == 5
        y = 100*(y-min(y))/min(y);
    end
    
    % optionally force weight at day 0 to 0
    forceDZeroWeight = false;
    % True: set weight to 0 at baseline;
    % False: don't set weight to 0 at baseline
    if forceDZeroWeight == true;
        y(1) = 0;  % force y(x==0) = 0;
    end

    % sort data
    [x, sortI] = sort(x);
    y = y(sortI);
    
    % assign final data to data structure
    pData(i).x = x;
    pData(i).y = y;
    pData(i).baselineWeight =  baselineWeight(i);   
    
    %  plot input data
    if ~isnan(pData(i).condition)
        subplot( p(2), p(1), pData(i).condition); hold on;
        plot (x, y);
         title([codingLegend.condition(pData(i).condition).label ': input data']);
    end
    
end

counter = 1;

% interpolate data (1 day intervals)
for i = 1:length(uMice)
    
    if max([pData(i).x]) >= max(xx)
        try;
            
            % determine interpolation range
            pData(i).xx = [min(xx):interpInterval:max(xx)];
            
            % interpolate data
            pData(i).yy = interp1([ pData(i).x], [ pData(i).y], [ pData(i).xx],  'pchip');
            
            % find Area under curve
            pData(i).AUC = trapz([pData(i).xx],[pData(i).yy]);

            % plot growth curves for each mouse
            figure(1);
            subplot( p(2), p(1), pData(i).condition); hold on;
            plot ([pData(i).xx], [pData(i).yy],'color', [0,0,0]+0.5);
            title(codingLegend.condition(pData(i).condition).label);
            
            % find rate of growth
            pData(i).xxDelta = pData(i).xx(2:end);
            pData(i).yyDelta = diff([pData(i).yy]);
            
            % plot rate of growth for each mouse
            figure(2);
            subplot( p(2), p(1), pData(i).condition); hold on;
            plot ([pData(i).xxDelta], [pData(i).yyDelta],'color', [0,0,0]+0.5);
            title(codingLegend.condition(pData(i).condition).label);
            
            % some organizational due diligence. 
            yDelta(:, i,pData(i).condition) = [pData(i).yyDelta];
            xDelta(:, i,pData(i).condition) = [pData(i).xxDelta];
            yAUC(:, i,pData(i).condition) = [pData(i).yy];
            xAUC(:, i,pData(i).condition) = [pData(i).xx];
            maxY(i) = max([ pData(i).yy]);
            minY(i) = min([ pData(i).yy]);
            maxYD(i) = max([ pData(i).yyDelta]);
            minYD(i) = min([ pData(i).yyDelta]);
        catch e; disp(e); end      
        
    else; 
        % specify which samples were omitted from analysis
        disp(['Sample ' num2str(uMice(i)) ' omitted due to insufficient data']);
    end
end

%% compute and plot average curves
for i = 1:length(unique([data.condition]))
    
    % retrive data from larger array
    yyRaw = yAUC(:,: ,i);
    yyDraw = yDelta(:,: ,i);
    
    % for each sample/mouse, consolidate data into more manageable format
    yyInput = [];     counter = 1;   
    for j = 1:size(yAUC,2)
        if sum(yyRaw(:,j)) ~= 0
            yyInput(:,counter) = yyRaw(:,j);
            yyDinput(:,counter) = yyDraw(:,j);
            counter = counter + 1;
        end
    end
    
    %growth data statistics
    yyMean(:,i) = mean(yyInput,2);
    yyStd(:,i) = std(yyInput,[],2);
    yySEM(:,i) = yyStd(:,i)/sqrt(size(yyInput,2));
    yyCIwidth(:,i) = yySEM(:,i)*1.96;
    
    %growth differential data statistics
    yyDMean(:,i) = mean(yyDinput,2);
    yyDStd(:,i) = std(yyDinput,[],2);
    yyDSEM(:,i) = yyDStd(:,i)/sqrt(size(yyDinput,2));
    yyDCIwidth(:,i) = yyDSEM(:,i)*1.96;
    
    % plot average growth curves
    figure(1);
    subplot( p(2), p(1), i); hold on;
    plot (xx, yyMean(:,i) ,'r', 'LineWidth',1.5);
    title(codingLegend.condition(i).label);
    xlabel('Time'); 
    
    % set y-axis labels
    if weightScale == 4||weightScale == 5
        ylabel('Growth (%)');
    elseif weightScale == 1||weightScale == 2||weightScale == 3
        ylabel('Growth');
    end
    
    % define y axis limits
    if min(minY) == 0
        ylim([0 max(maxY)+(0.1*max(maxY))]);
    else
        ylim([min(minY)-(0.1*max(maxY))  max(maxY)+(0.1*max(maxY))]);
    end
    
     % plot average growth rate curves
    figure(2);
    subplot( p(2), p(1), i); hold on;
    plot (xx(2:end), yyDMean(:,i) ,'r', 'LineWidth',1.5);
    title(codingLegend.condition(i).label);
    xlabel('Time'); 
    
    % set y-axis labels
    if  weightScale == 4||weightScale == 5
        ylabel('Rate of Growth (% per unit time)');
    elseif weightScale == 1||weightScale == 2||weightScale == 3
        ylabel('Rate of Growth (per unit time)');
    end   
    
    % define y axis limits
    ylim([min(minYD)-(0.1*max(maxYD))  max(maxYD)+(0.1*max(maxYD))]);
    
    % plot no growth reference line
    hline(0,'k:')
end

%% select subset of data to visualize side by side.

% specify available options
optns = {codingLegend.condition.label};
optns{end+1} = 'all of the above';
optns{end+1} = 'continue with analysis...';
fin = 0;
choice = [];

% allow user input to select which curves to compare
if length(optns)>4
    while fin == 0
        choice = [choice menu('which subset of conditions would you like to compare directly?',  optns)];
        disp(optns{choice(end)})
        if choice(end) == length(optns) | choice(end) == length(optns)-1
            fin = 1;
        else
        end
    end
else
    choice = length(optns)-1; % if there are only 2 options, select both for comparison
end

% some mental gymnastics
if choice(end) == length(optns)
    choice = choice(1:end-1);
elseif choice(end) == length(optns)-1;
    choice = 1:length(unique([data.condition]));
elseif choice(end) == length(optns) && length(choice) == 1
    choice = 1:length(unique([data.condition]));
end

% specify transparency coeffiicent that wil be used for 95% CI curves
if length(choice) < 3
    transparency = 0.35;
else
    transparency = 0.20;
end

% plot subset of curves specified in above menu
subset1 = []; subset2 = [];
counter = 1;
for i = 1:length(unique([data.condition]))
    if sum(choice == i)>0
        % plot subset of average growth curves
        figure(3);
        h1{counter} = plot (xx, yyMean(:,i)); hold on;
        color{counter} = get(h1{counter}, 'Color');
        ciplot(yyMean(:,i)-yyCIwidth(:,i),yyMean(:,i)+yyCIwidth(:,i), xx, color{counter}, transparency);
        plotLegend{counter} = codingLegend.condition(i).label;
        subset1 = [subset1 h1{counter}(1)];
        
        %plot subset of average growth rate curves
        figure(4);
        h1{counter} = plot (xx(2:end), yyDMean(:,i)); hold on;
        ciplot(yyDMean(:,i)-yyDCIwidth(:,i),yyDMean(:,i)+yyDCIwidth(:,i), xx(2:end), color{counter}, transparency);
        subset2 = [subset2 h1{counter}(1)];
        
        counter = counter + 1;
    end
end

% assign legend, title, x and y axis labels. 
figure(3); legend(subset1, plotLegend);  title('Growth Curves (Average)');  xlabel('Time');

if  weightScale == 4||weightScale == 5
    ylabel('Growth (%)');
elseif weightScale == 1||weightScale == 2||weightScale == 3
    ylabel('Growth');
end

figure(4); legend(subset2, plotLegend); title('Growth Rate Curves (Average)'); xlabel('Time'); 

if  weightScale == 4||weightScale == 5
    ylabel('Rate of Growth (% per unit time)');
elseif weightScale == 1||weightScale == 2||weightScale == 3
    ylabel('Rate of Growth (per unit time)');
end

% assign no growth reference line
hline(0,'k:')

% determine range of AUC values computed, to be used later for figure axis specification
minAUC = min([pData.AUC]);
maxAUC = max([pData.AUC]);

%% results and distribution of results

figure(5);
for i = 1:length(unique([data.condition]))
    
    % assign current condition to results data structure.
    results(i).condition = (codingLegend.condition(i).label);
    
    % assign AUC data to temporary variable for easier manipulation
    try;
        tempData1 = [pData([pData.condition]==i).AUC]';
    catch; error('check if range of dates used in analysis is too wide'); end
    IDmice = [pData([pData.condition]==i).ID]';
    
    % assign baseline weight and age to results
    try;
        DOB = [pData([pData.condition]==i).birthDate]';
        SOD = [pData([pData.condition]==i).startDate]';
        baselineWeight = [pData([pData.condition]==i).baselineWeight]';
        age = SOD-DOB;
        results(i).age = age([tempData1]>0);
        results(i).baselineWeight = baselineWeight([tempData1]>0);
    catch e; disp(e); end
    
    % assign sample ID, AUC, mean AUC and median AUC to results
    results(i).ID = IDmice([tempData1]>0);
    results(i).AUC = tempData1([tempData1]>0);
    results(i).meanAUC = mean([results(i).AUC]);
    results(i).medianAUC = median([results(i).AUC]);
    
    % plot histogram of computed AUC values  
    yLabels{i} =  results(i).condition;
    subplot(p(2), p(1), i);
    hist( [results(i).AUC] );
    xlim([minAUC maxAUC])
    title(codingLegend.condition(i).label);
    xlabel('AUC'); ylabel('Count');
end

% plot average AUC in bar plot
figure (6);
barh([results.meanAUC ],  'FaceColor', [0.75 0.75 0.75]); hold on;
vline(mean([results.meanAUC ],'omitnan'),'r'); % plot overall mean cv
set(gca,'yticklabel', yLabels);
xlabel('AUC');
title('Growth Curve AUC');


%% save results

% save results to excel file only if saveResults == TRUE
if saveResults
    for i = 1:length(results)
        results(i).condition(results(i).condition == ' ') = [];
        for j = 1:length([results(i).ID])
            try; ageRes(j).(results(i).condition) = results(i).age(j); catch e; disp(e);end
            try; weightRes(j).(results(i).condition) = results(i).baselineWeight(j); catch e; disp(e);end
            try;  AUCres(j).(results(i).condition) = results(i).AUC(j); catch e; disp(e);end
            try; IDres(j).(results(i).condition) = results(i).ID(j); catch e; disp(e); end
        end
    end
    try;
        % get current time
        now = datestr(datetime('now'));
        now = strrep(now, ':','-');
        
        % assign file name for results
        file = strrep(file, '.xlsx',' ');
        file = [file 'RESULTS ' now '.xlsx'];
        
        % save results to excel file
        writetable(struct2table(IDres), file, 'Sheet' ,['ID']);
        writetable(struct2table(AUCres), file, 'Sheet' ,['AUC']);
        writetable(struct2table(ageRes), file, 'Sheet' ,['age']);
        writetable(struct2table(weightRes), file, 'Sheet' ,['baseline Weight']);
        
    catch;
        error(' ')
    end
    % remove first 3 blank sheets in new excel file (Sheet 1, Sheet 2, Sheet 3)
    RemoveSheet123(file);
end

% inform user that analysis is complete!
display('Analysis Complete!');



%% functions
function [data] = date2time(data, dateVar)

for i = 1:length(dateVar)
    for j = 1:length(data)
        try;
            data(j).(dateVar{i}) = datenum(data(j).(dateVar{i}));
        catch; data(j).(dateVar{i}) = nan(); end
    end
end

end

function [data, codingLegend] = cat2num(data, categoricalVar);
for i = 1:length(categoricalVar)
    
    % define coding scheme
    a= categorical(cellstr({data.(categoricalVar{i})}));
    [b,GN,~]= grp2idx(a);
    
    % create coding legend
    for j = 1:length(GN)
        codingLegend.(categoricalVar{i})(j).code = j;
        codingLegend.(categoricalVar{i})(j).label = GN{j};
    end
    
    % assign codes to dataset
    for j = 1:length(data)
        data(j).(categoricalVar{i}) = b(j);
    end
    
end
end

function data = importData(file, sheet)

% import spreadsheet
[num,txt,raw]=xlsread(file,sheet);

% extract data and store in data structure
for i = 2:size(raw,1);
    for j = 1:size(raw,2)
        if ~isnan(raw{1,j}); data(i-1).(raw{1,j}) = raw{i,j}; end
    end; end
end


function [p,n]=numSubplots(n)
% function [p,n]=numSubplots(n)
%
% Purpose
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots.
%
% Inputs
% n - the desired number of subplots.
%
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
%
%
% Example: neatly lay out 13 sub-plots
% >> p=numSubplots(13)
% p =
%     3   5
% for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end
%
%
% Rob Campbell - January 2010


while isprime(n) & n>4,
    n=n+1;
end
p=factor(n);
if length(p)==1
    p=[1,p];
    return
end
while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end
    p=sort(p);
end
%Reformat if the column/row ratio is too large: we want a roughly
%square design
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end
end

function [fitresult, gof] = fitHill(x, y)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 12-Jan-2019 12:17:52


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'b1*(x^b2)/((b3^b2)+(x^b2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.Upper = [inf inf inf];
opts.StartPoint = [0.957506835434298 25 30];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel x
% ylabel y
% grid on
end

function hhh=vline(x,in1,in2)
% Credit: Brandon Kuczenski
% https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline
%
% function h=vline(x, linetype, label)
%
% Draws a vertical line on the current axes at the location specified by 'x'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = vline(42,'g','The Answer')
%
% returns a handle to a green vertical line on the current axes at x=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% vline also supports vector inputs to draw multiple lines at once.  For example,
%
% vline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
%
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001
if length(x)>1  % vector input
    for I=1:length(x)
        switch nargin
            case 1
                linetype='r:';
                label='';
            case 2
                if ~iscell(in1)
                    in1={in1};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                label='';
            case 3
                if ~iscell(in1)
                    in1={in1};
                end
                if ~iscell(in2)
                    in2={in2};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                if I>length(in2)
                    label=in2{end};
                else
                    label=in2{I};
                end
        end
        h(I)=vline(x(I),linetype,label);
    end
else
    switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            linetype=in1;
            label='';
        case 3
            linetype=in1;
            label=in2;
    end
    
    
    
    g=ishold(gca);
    hold on
    y=get(gca,'ylim');
    h=plot([x x],y,linetype);
    if length(label)
        xx=get(gca,'xlim');
        xrange=xx(2)-xx(1);
        xunit=(x-xx(1))/xrange;
        if xunit<0.8
            text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        else
            text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        end
    end
    if g==0
        hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end % else
if nargout
    hhh=h;
end

end


function ciplot(lower,upper,x,colour, transparency);

% ciplot(lower,upper)
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.
% Raymond Reynolds 24/11/06
if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end
if nargin<4
    colour='b';
end
if nargin<3
    x=1:length(lower);
end
% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
    x=x'; end
if find(size(lower)==(max(size(lower))))<2
    lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
    upper=upper'; end
h1 = fill([x fliplr(x)],[upper fliplr(lower)],colour)
set(h1,'EdgeColor','none')
set(h1, 'facealpha', transparency);
end


function hhh=hline(y,in1,in2)
% function h=hline(y, linetype, label)
%
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% hline also supports vector inputs to draw multiple lines at once.  For example,
%
% hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
%
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001
if length(y)>1  % vector input
    for I=1:length(y)
        switch nargin
            case 1
                linetype='r:';
                label='';
            case 2
                if ~iscell(in1)
                    in1={in1};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                label='';
            case 3
                if ~iscell(in1)
                    in1={in1};
                end
                if ~iscell(in2)
                    in2={in2};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                if I>length(in2)
                    label=in2{end};
                else
                    label=in2{I};
                end
        end
        h(I)=hline(y(I),linetype,label);
    end
else
    switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            linetype=in1;
            label='';
        case 3
            linetype=in1;
            label=in2;
    end
    
    
    
    g=ishold(gca);
    hold on
    x=get(gca,'xlim');
    h=plot(x,[y y],linetype);
    if ~isempty(label)
        yy=get(gca,'ylim');
        yrange=yy(2)-yy(1);
        yunit=(y-yy(1))/yrange;
        if yunit<0.2
            text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
        else
            text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
        end
    end
    if g==0
        hold off
    end
    set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
end % else
if nargout
    hhh=h;
end

end

function RemoveSheet123(excelFileName,sheetName)
% RemoveSheet123 - removes the sheets that are automatically added to excel
% file.
% When Matlab writes data to a new Excel file, the Excel software
% automatically creates 3 sheets (the names are depended on the user
% languade). This appears even if the user has defined the sheet name to be
% added.
%
% Usage:
% RemoveSheet123(excelFileName) - remove "sheet1", "sheet2","sheet3" from
% the excel file. excelFileName is a string of the Excel file name.
% RemoveSheet123(excelFileName,sheetName) - enables the user to enter the
% sheet name when the language is other than English.
% sheetName is the default sheet name, without the number.
%
%
%                       Written by Noam Greenboim
%                       www.perigee.co.il
%
%% check input arguments
if nargin < 1 || isempty(excelFileName)
    error('Filename must be specified.');
end
if ~ischar(excelFileName)
    error('Filename must be a string.');
end
try
    excelFileName = validpath(excelFileName);
catch
    error('File not found.');
end
if nargin < 2
    sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, HE: ?????? , etc. (Lang. dependent)
else
    if ~ischar(sheetName)
        error('Default sheet name must be a string.');
    end
end
%%
% Open Excel file.
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(excelFileName); % Full path is necessary!
% Delete sheets.
try
    % Throws an error if the sheets do not exist.
    objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
    fprintf('\nsheet #1 - deleted.')
    objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
    fprintf('\nsheet #2 - deleted.')
    objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
    fprintf('\nsheet #3 - deleted.\n')
catch
    fprintf('\n')
    O=objExcel.ActiveWorkbook.Worksheets.get;
    if O.Count==1
        error('Can''t delete the last sheet. Excel file must containt at least one sheet.')
    else
        warning('Problem occured. Check excel file.');
    end
end
% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;
end
function filenameOut = validpath(filename)
% VALIDPATH builds a full path from a partial path specification
%   FILENAME = VALIDPATH(FILENAME) returns a string vector containing full
%   path to a file. FILENAME is string vector containing a partial path
%   ending in a file or directory name. May contain ..\  or ../ or \\. The
%   current directory (pwd) is prepended to create a full path if
%   necessary. On UNIX, when the path starts with a tilde, '~', then the
%   current directory is not prepended.
%
%   See also XLSREAD, XLSWRITE, XLSFINFO.

%   Copyright 1984-2012 The MathWorks, Inc.

%First check for wild cards, since that is not supported.
if strfind(filename, '*') > 0
    error(message('MATLAB:xlsread:Wildcard', filename));
end

% break partial path in to file path parts.
[Directory, file, ext] = fileparts(filename);
if ~isempty(ext)
    filenameOut = getFullName(filename);
else
    extIn = matlab.io.internal.xlsreadSupportedExtensions;
    for i=1:length(extIn)
        try                                                                %#ok<TRYNC>
            filenameOut = getFullName(fullfile(Directory, [file, extIn{i}]));
            return;
        end
    end
    error(message('MATLAB:xlsread:FileDoesNotExist', filename));
end
end
function absolutepath=abspath(partialpath)

% parse partial path into path parts
[pathname, filename, ext] = fileparts(partialpath);
% no path qualification is present in partial path; assume parent is pwd, except
% when path string starts with '~' or is identical to '~'.
if isempty(pathname) && strncmp('~', partialpath, 1)
    Directory = pwd;
elseif isempty(regexp(partialpath,'(.:|\\\\)', 'once')) && ...
        ~strncmp('/', partialpath, 1) && ...
        ~strncmp('~', partialpath, 1);
    % path did not start with any of drive name, UNC path or '~'.
    Directory = [pwd,filesep,pathname];
else
    % path content present in partial path; assume relative to current directory,
    % or absolute.
    Directory = pathname;
end

% construct absolute filename
absolutepath = fullfile(Directory,[filename,ext]);
end
function filename = getFullName(filename)
FileOnPath = which(filename);
if isempty(FileOnPath)
    % construct full path to source file
    filename = abspath(filename);
    if isempty(dir(filename)) && ~isdir(filename)
        % file does not exist. Terminate importation of file.
        error(message('MATLAB:xlsread:FileDoesNotExist', filename));
    end
else
    filename = FileOnPath;
end
end