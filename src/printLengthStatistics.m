function printLengthStatistics(tValues,acc2)

tLengths = tValues(:,2) - tValues(:,1);

accTs = abs(tValues(acc2,3));
rejTs = abs(tValues(~acc2,3));


lengths = {
    tLengths
    abs(tValues(:,3))
    accTs
    rejTs
};

names = {
    'Line segments'
    'Proposed relocations'
    'Accepted relocations'
    'Rejected relocations'
};

stats = cell2mat(cellfun(@(x) ...
    [mean(x), median(x), min(x), max(x)], ...
    lengths, 'UniformOutput', false));

T = array2table(stats, ...
    'VariableNames', {'Mean', 'Median', 'Min', 'Max'}, ...
    'RowNames', names);

disp('Length statistics, Node relocations:')
disp(newline)
disp(T)

end
