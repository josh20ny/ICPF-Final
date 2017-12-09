%Runs the scripts that plots the data into multiple images 
for n = 10:2:20
    run(['GraphAt' num2str(n) 'Days.m']);
end

%Runs the script that plots the graph with error bars

%Runs the script that plots the Gompertzian fit along side the observed graph
run('mathematicalModel.m');