# ICPF-Final

This is the Final project for COE 301 Introduction to Computer Programming Fall 2017 by Josh, Darcey and Rhythem.

---

# Note in order to run properly with main.m

In order to make sure the code runs properly, make sure that whichever folder you are in, when you run **main.m**, is organized as follows:

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/Organization%20Schematics.png)

You can pull all of the source codes from each of the files in **src** and the **cells.mat** from the **data** file.

---

# Main

The entire project is run through this script, main.m:

```MATLAB

%Runs the scripts that plots the data into multiple images 
for n = 10:2:20
    run(['GraphAt' num2str(n) 'Days.m']);
end

%Runs the script that plots the graph with error bars
run('errorInCellCount.m');

%Runs the script that plots the Gompertzian fit along side the observed graph
run('mathematicalModel.m');
```

# Data Reduction and Visualization (Darcey Long)

Visualization of the Tumor Growth over 20 days

GraphAt10Days.m (the code is similar for the other day's scripts)

```MATLAB
clear all;
load('cells.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE SPECIFICATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrow = round(sqrt(length(cells(1,1,:,:)))); % number of subplots to be created in the y direction
ncol = nrow;                             % number of subplots to be created in the x direction

% main plot specifications
mainPlotMarginTop = 0.06;       % top margin of the main axes in figure
mainPlotMarginBottom = 0.12;    % bottom margin of the main axes in figure
mainPlotMarginLeft = 0.08;      % bottom margin of the main axes in figure
mainPlotMarginRight = 0.1;      % right margin of the main axes in figure
mainPlotPositionX = 0.05;       % the x coordinate of the bottom left corner of the main axes in figure
mainPlotPositionY = 0.08;       % the y coordinate of the bottom left corner of the main axes in figure
mainPlotWidth = 1 - mainPlotMarginRight - mainPlotPositionX; % the width of the main axes in figure
mainPlotHeight = 1 - mainPlotMarginTop - mainPlotPositionY; % the height of the main axes in figure
mainPlotTitleFontSize = 12;     % The font size for the main plot labels and title
mainPlotAxisFontSize = 12;      % The font size for the main plot labels and title

% subplot properties
subPlotFontSize = 10;     % the font size for subplots
subplotInterspace = 0.03; % space between subplots
subplotWidth = (1-mainPlotMarginLeft-mainPlotMarginRight-nrow*subplotInterspace)/ncol;   % The width of subplots
subplotHeight = (1-mainPlotMarginTop-mainPlotMarginBottom-ncol*subplotInterspace)/nrow ; % The height of subplots

% specifications of the color bar to the figure
colorbarFontSize = 13;                                           % the font size of the color bar
colorbarWidth = 0.03;                                            % the width of the color bar
colorbarPositionY = mainPlotMarginBottom;                        % the y-position of the color bar
colorbarPositionX = 1 - mainPlotMarginRight;                     % the x-position of the color bar
colorbarHeight = nrow*subplotHeight+(nrow-1)*subplotInterspace;  % the height of the color bar
colorLimits = [0,max(max(max(cells(:,:,:))))];                    % the color bar limits, the dataset contains numbers between 0 and some large number.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST CREATE A FIGURE HANDLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figHandle = figure();                           % create a new figure
figHandle.Visible = 'on';                       % set the visibility of figure in MATLAB
figHandle.Position = [0, 0, 900, 1350];         % set the position and size of the figure
figHandle.Color = [0.9400 0.9400 0.9400];       % set the background color of the figure to MATLAB default. You could try other options, such as 'none' or color names 'red',...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD THE MAIN AXES TO THE FIGURE: 
% The main axes is needed in order to add
% the x and y labels and the color bar
% for the entire figure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mainPlot = axes();              % create a set of axes in this figure. This main axes is needed in order to add the x and y labels and the color bar for the entire figure.
mainPlot.Color = 'none';        % set color to none
mainPlot.FontSize = 11;         % set the main plot font size
mainPlot.Position = [ mainPlotPositionX mainPlotPositionY mainPlotWidth mainPlotHeight ]; % set position of the axes
mainPlot.XLim = [0 1];          % set X axis limits
mainPlot.YLim = [0 1];          % set Y axis limits
mainPlot.XLabel.String = 'Voxel Number in X Direction'; % set X axis label
mainPlot.YLabel.String = 'Voxel Number in Y Direction'; % set Y axis label
mainPlot.XTick = [];            % remove the x-axis tick marks
mainPlot.YTick = [];            % remove the y-axis tick marks
mainPlot.XAxis.Visible = 'off'; % hide the x-axis line, because we only want to keep the x-axis label
mainPlot.YAxis.Visible = 'off'; % hide the y-axis line, because we only want to keep the y-axis label
mainPlot.XLabel.Visible = 'on'; % make the x-axis label visible, while the x-axis line itself, has been turned off;
mainPlot.YLabel.Visible = 'on'; % make the y-axis label visible, while the y-axis line itself, has been turned off;
mainPlot.Title.String = ['Time = 10.0 Days. Brain MRI slices along Z-direction, Rat W09. No radiation treatment']; % set the title of the figure
mainPlot.XLabel.FontSize = mainPlotAxisFontSize; % set the font size for the x-axis in mainPlot
mainPlot.YLabel.FontSize = mainPlotAxisFontSize; % set the font size for the y-axis in mainPlot
mainPlot.Title.FontSize = mainPlotTitleFontSize; % set the font size for the title in mainPlot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD COLORBAR TO THE FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now add the color bar to the figure
axes(mainPlot);                       % focus on the mainPlot axes, because this is where we want to add the colorbar
caxis(colorLimits);                   % set the colorbar limits
cbar = colorbar;                      % create the color bar!
ylabel(cbar,'Number of Tumor Cells'); % now add the color bar label to it
cbar.Position = [ colorbarPositionX ... Now reset the position for the colorbar, to bring it to the rightmost part of the plot
                  colorbarPositionY ...
                  colorbarWidth ...
                  colorbarHeight ...
                ];
cbar.FontSize = colorbarFontSize;     % set the fontsize of the colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD THE SUBPLOTS TO THE FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sliceNumber = 0;
for irow = nrow:-1:1
    for icol = 1:ncol
        sliceNumber = sliceNumber + 1;
        subPlot = axes( 'position', [ ... set the position of the axes for each subplot
                                      (icol-1)*(subplotInterspace+subplotWidth) + mainPlotMarginLeft ...
                                      (irow-1)*(subplotInterspace+subplotHeight) + mainPlotMarginBottom ...
                                      subplotWidth ...
                                      subplotHeight ...
                                    ] ...
                      );
        imagesc(cells(:,:,sliceNumber),colorLimits);
        hold on;
        if icol~= 1
            subPlot.YTickLabels = [];
        end
        if irow ~= 1
            subPlot.XTickLabels = [];
        end
     subPlot.Title.String = ['z = ', num2str(sliceNumber)];
     subPlot.CLim = colorLimits;
     
        % THE REST OF THIS CODE SHOULD BE ADDED BY YOU:
        % Here you put the image slice that you want to plot (call imagesc() with the appropriate input data slice)
        % then hold on to this figure.
        % Now call the relevant properties of subplot to modify it properties, for example,
        % to ensure that whenever the subplot is NEITHER in the first columns of plots
        % NOR in the bottom row of subplots in the figure, then its XTickLabels and YTickLabels
        % are set to null value []. For this, you will need to write an if-block.
        % Also add the title for each subplot here using the subplot property subPlot.Title.String
        % Finally set the color limits of the subplot to the proper range by subPlot.CLim = colorLimits
        % where the variable colorLimits is defined in the above lines when we created the colorbar.
        % END OF YOUR CODE ADDITIONS
    end
end

saveas(gcf,'PlotAtDay10.png');        % save the figure
```

**Day 10**

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/results/Data%20Visualization%20(Part%201)/PlotAtDay10.png)

**Day 12**

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/results/Data%20Visualization%20(Part%201)/PlotAtDay12.png)

**Day 14**

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/results/Data%20Visualization%20(Part%201)/PlotAtDay14.png)

**Day 16**

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/results/Data%20Visualization%20(Part%201)/PlotAtDay16.png)

**Day 18**

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/results/Data%20Visualization%20(Part%201)/PlotAtDay18.png)

**Day 20**

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/results/Data%20Visualization%20(Part%201)/PlotAtDay20.png)

---

# Obtaining the Error in Tumor Cell Count (Joshua Montoya)

errorInCellCount.m

```MATLAB
% Loading the experiment Data
load('cells.mat');

%initilizing an array for the time, total Cells, and error in the Cell
%Growth
totalCells = {zeros(1 + length(cells(1,1,1,:)),1),zeros(1 + length(cells(1,1,1,:)),1)};
errorArray = zeros( (length(cells(1, 1, 1, :)) + 1), 1);

%Sets the days, sums the total cells, and sums the error cells
for day = 1:length(totalCells{1})
    if day == 1
        totalCells{1}(day) = 0;
        totalCells{2}(day) = 0;
        errorArray(day) = 0;
    else
        totalCells{1}(day) = 10 + 2*(day-2);
        for slicenumber = 1:length(cells(1, 1, :, 1))
            BW = imbinarize(cells(:, :, slicenumber, day-1));
            [B, L] = bwboundaries(BW, 'noholes');
            if ~isempty(B)
                for sizeB = 1:length(B)
                    Array = B{sizeB};
                    for q = 1:(length(Array) - 1)
                        cellError = cellError + sum(sum(cells(Array(q, 1), Array(q, 2), slicenumber, day-1)));
                    end
                end
            end
            totalCells{2}(day) = totalCells{2}(day) + sum(sum(cells(:,:,slicenumber,(day-1))));
        end
        errorArray(day) = cellError;
    end
    cellError = 0;
end

% displays the data to be plotted
disp(totalCells{1});
disp(totalCells{2});
disp(errorArray);

% Plots the data and error
figure();
errorbar(totalCells{1},totalCells{2},errorArray, '.-b', 'linewidth' ... 
    , 4, 'MarkerSize', 32);

title('Rat Brain Tumor Cell Growth');
xlabel('Time [Days]');
ylabel('Experimental Tumor Cell Count');
legend('Experimental Data', 'Location', 'Northwest')
box on;

saveas(gcf, 'Cell Error.png');
```

# Graph

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/results/Error%20(Part%202)/Cell%20Error.png)

---

# Mathematical Modeling (Rhythem Sharma)

This is the code used to develop the best fit line as shown in the graph below.

mathematicalModel.m

```MATLAB
load('cells.mat');
global timeArray dataVals;
dataVals = zeros(8, 1);
dataVals(1) = 100000;
timeArray = [0, 10 12 14 16 18 20 22];
totalCellError = 0;
errArray = zeros(8, 1);
errArray(1) = 0;
params = [10, 0.1, 1]; % initial values
tint = 0:.1:25;

GompGrowth = @(t, N, lamda, c) N*exp(lamda*(1-exp(-c*t)));

%includes finding the number of cells along the side
for k = 1:7
    for n = 1:16
        BI = imbinarize(cells(:, :, n, k));
        [B, L] = bwboundaries(BI, 'noholes');
        if(~isempty(B))
            for m = 1:length(B)
                Array = B{m};
                disp([num2str(k) ' ' num2str(n)]);
                disp(Array)
                for x = 1:length(Array) - 1
                    disp(totalCellError)
                    totalCellError = totalCellError + cells(Array(x, 1), Array(x, 2), n, k);
                end
            end
        end
        dataVals(k + 1) = dataVals(k + 1) + sum(sum(cells(:, :, n, k)));
    end
    errArray(k + 1) = totalCellError/2;
    totalCellError = 0;
end

errorbar(timeArray, dataVals, errArray, 'b.-', 'linewidth', 3, 'MarkerSize', 32, 'DisplayName', 'Experimental Data');
hold on;

%summing the valus
optimizedParamVals = fminsearch(@sumLogs, params);
plot(tint, GompGrowth(tint, dataVals(1), optimizedParamVals(1), optimizedParamVals(2)), 'linewidth', 4, 'color', 'r', 'DisplayName', 'Gomertizan Fit');

%formatting the plotted figure
xlabel('Time [Days]');
ylabel('Tumor Cell Count');
title('Gompertzian Fit to Rat''s Brain Tumor Growth', 'fontsize', 15);
legend('show', 'location', 'northwest');
saveas(gcf, 'Gompertizan Fit plot.png');

%writing all the optimized values to a file
fileParams = fopen('fileWithParams.txt', 'w');
fprintf(fileParams, ['Lamda: ', num2str(optimizedParamVals(1)), ...
    ', C-Value: ', num2str(optimizedParamVals(2)), ...
    ', Sigma: ', num2str(optimizedParamVals(3))]);
fclose(fileParams);

%using a function to find the probability density function value
function output = sumLogs(params)
    global dataVals timeArray;
    
    %physical model
    GompGrowth = @(t, N, lamda, c) N*exp(lamda*(1-exp(-c*t)));

    %probability density model
    probDensity = @(N_Obs, sigma, GompGrowthVal) log((1/(N_Obs*sigma*sqrt(2*pi))) * ...
        exp((-(log(N_Obs) - log(GompGrowthVal))^2)/(2*sigma^2)));
    
    output = 0;
    for n = 1:length(dataVals)
        GompGrowthVal = GompGrowth(timeArray(n), dataVals(1), params(1), params(2));
        output = output + probDensity(dataVals(n), params(3), GompGrowthVal);
    end
    output = -1 * output;
end
```

# Graph

![Image](https://github.com/josh20ny/ICPF-Final/blob/master/results/Gompertzian%20Fit%20(Part%203)/Gompertizan%20Fit%20plot.png)

The 3 specific values that were required in order for this best fit line to be modeled can be found in the **fileWithParams.txt** file. In there are these three values:

`Lamda: 8.0729, C-Value: 0.1102, Sigma: 0.17728`
