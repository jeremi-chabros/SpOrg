clearvars;
load lookup2

% top 4 rows: 12 13 14 15 16 17 21 22 23 24 25 26 27 28 31 32 33 34 35 36
% 37 38 41 42 43 44 45 46 47 48
% col1 [12 13 14 16 17];
% col8 [82 83 84 85 86 87];
% col12 [12,13,14,16,17,21,22,23,24,25,26,27,28];
% col78 [71 72 73 74 75 76 77 78 82 83 84 85 86 87];

files(1).stim = [82 83 84 85 86 87];
files(2).stim = [12,13,14,16,17,21,22,23,24,25,26,27,28];
files(3).stim = [];%
files(4).stim = [82 83 84 85 86 87];
files(5).stim = [12,13,14,16,17,18,21,22,23,24,25,26,27,28];
files(6).stim = [71 72 73 74 75 76 77 78 82 83 84 85 86 87];
files(7).stim = [12,13,14,16,17,21,22,23,24,25,26,27,28];
files(8).stim = [];%
files(9).stim = [];%
files(10).stim = [12,13,14,16,17,21,22,23,24,25,26,27,28];
files(11).stim = [12 13 14 16 17];
files(12).stim = [12 13 14 16 17];
files(13).stim = [82 83 84 85 86 87];
files(14).stim = [12 13 14 16 17];
files(15).stim = [82 83 84 85 86 87];
files(16).stim = [82 83 84 85 86 87];
files(17).stim = [12 13 14 16 17];
files(18).stim = [];%
files(19).stim = [12 13 14 16 17];
files(20).stim = [];%
files(21).stim = [];%
files(22).stim = [82 83 84 85 86 87];
files(23).stim = [12,13,14,16,17,21,22,23,24,25,26,27,28];
files(24).stim = [];%
files(25).stim = [];%
files(26).stim = [12 13 14 16 17];
files(27).stim = [];%
files(28).stim = [12,13,14,16,17,21,22,23,24,25,26,27,28];

save('lookup2.mat', 'files');