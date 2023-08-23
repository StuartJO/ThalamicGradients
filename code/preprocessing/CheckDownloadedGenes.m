function CheckDownloadedGenes()

% This is quite frankly a terrible way to check which genes were actually
% downloaded, but it was how I did it as a first pass which became the only
% pass. I have kept it as is because the following returns the genes in a
% specific order and that order is used throughout

GeneDir = './data/AHBA_wholebrain/Genes';

GeneDir_dir = dir(GeneDir);

GeneDir_dir(1:2) = [];

GenesKept = zeros(length(GeneDir_dir),1);

for i = 1:length(GeneDir_dir)
    GenesKept(i) = str2double(gene_dir(i).name);
end

save('./data/preprocessed/GenesKept.mat','GenesKept')
dlmwrite('./data/preprocessed/GenesKept.txt',GenesKept','Precision',9)