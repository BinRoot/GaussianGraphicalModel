load('../data/entrez.txt');
res1_indices = importdata('res1_indices.mat');
genePairs = arrayfun(@(x) entrez(x), res1_indices);
disp(genePairs)
