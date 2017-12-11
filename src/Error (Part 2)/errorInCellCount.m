cd ../../data;

load('cells.mat');

cd ../src/'error (Part 2)';

totalCells = {zeros(1 + length(cells(1,1,1,:)),1),zeros(1 + length(cells(1,1,1,:)),1)};
errorArray = zeros( (length(cells(1, 1, 1, :)) + 1), 1);

for day = 1:length(totalCells{1})
    if day == 1
        totalCells{1}(day) = day;
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
        errorArray(day) = cellError/2;
    end
    cellError = 0;
end

disp(totalCells{1});
disp(totalCells{2});
disp(errorArray);

errorbar(totalCells{1},totalCells{2},errorArray, 'b', 'linewidth' ... 
    , 4, 'MarkerSize', 32);

title('Rat Brain Tumor Cell Growth');
xlabel('Time [Days]');
ylabel('Experimental Tumor Cell Count');
legend('Experimental Data', 'Location', 'Northwest')
box on;

cd ../../results/'Error (Part 2)'

saveas(gcf, 'Cell Error.png');