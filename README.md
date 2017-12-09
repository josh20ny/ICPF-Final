# ICPF-Final

This is the Final project for COE 301 Introduction to Computer Programming Fall 2017 by Josh, Darcey and Rhythem.

---

# Data Reduction and Visualization

Visualization of the Tumor Growth over 20 days

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

# Obtaining the Error in Tumor Cell Count

---

# Mathematical Modeling

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

# main.m

This script runs the entire project
