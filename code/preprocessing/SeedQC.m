function SeedQC()

% This function performs quality control (QC) on seed data, filters seeds,
% and processes tractography data.

type = 'fsl';
mm = '1.75';

% Set the threshold for the average correlation between subject's thalamic seed vectors
SubThr = .7;

% Set the threshold for the proportion of seeds across good QC'd subjects
% that exist in individual thalamic masks
SeedThr = .85;

% Load the list of unrelated subjects
SUB = dlmread('UnrelatedSubs.txt');

% Determine the number of seeds
Nseeds = length(dlmread(['./data/preprocessed/seeds_',mm,'mm_ind.txt']));
DATA = zeros(Nseeds, length(SUB)); % Initialize a matrix to store seed data

% Populate the DATA matrix with subject's seed data
for i = 1:length(SUB)
    DATAtemp = dlmread(['./data/tractography/SUBJECTS/', num2str(SUB(i)), '/masked_seed_ind_', mm, 'mm_', type, '_rows_', num2str(SUB(i)), '.txt']);
    DATA(DATAtemp, i) = 1;
end

% Load and filter gene-specific seed data
InGene = dlmread('gene_masked_seeds_1.75mm_ind.txt');
InGeneLogical = false(Nseeds, 1);
InGeneLogical(InGene) = true;
DATA(~InGeneLogical, :) = 0;

% imagesc(DATA)
% ylabel('Seed number')
% xlabel('Subject')
% set(gca,'FontSize',20)
% title('Yellow = in mask, blue = outside mask')

% Calculate correlation matrix
C = corr(DATA);

% figure
% imagesc(C)
% cc = colorbar;
% cc.Label.String = 'Correlation';
% set(gca,'FontSize',20)

% Calculate the mean correlation for each seed
mean_corr = mean(C);

% Store original seed data
DATA_ORIG = DATA;

% Determine valid seeds based on average correlation threshold
DATA_KEEP_IND = mean_corr >= SeedThr;
GOOD_SUB = SUB(DATA_KEEP_IND);

% Filter seed data and calculate proportions
DATA = DATA_ORIG(:, DATA_KEEP_IND);
SubValid = sum(DATA, 1) ./ Nseeds;
NSubValid = length(SubValid);
SeedValid = sum(DATA, 2) ./ length(GOOD_SUB);
ValidSeed = SeedValid > SubThr;

% Save filtered seed data and information
dlmwrite(['./data/preprocessing/seeds_', mm, 'mm_', type, '_seeds2use_', num2str(NSubValid), 'subs_valid_bin.txt'], double(ValidSeed)')
dlmwrite(['./data/preprocessing/seeds_', mm, 'mm_', type, '_seeds2use_', num2str(NSubValid), 'subs_valid_thr.txt'], SeedValid')

% Process tractography data for each good subject
for i = 1:length(GOOD_SUB)
    tracks = read_mrtrix_tracks(['./data/tractography/SUBJECTS/', num2str(GOOD_SUB(i)), '/', num2str(GOOD_SUB(i)), '_seeds_', mm, 'mm_fsl.tck']);
    tracts = tracks.data;
    valid_seeds = tracts(ValidSeed)';
    allGeneSeeds = tracts(InGeneLogical)';
    dlmwrite(['./data/tractography/SUBJECTS/', num2str(GOOD_SUB(i)), '/', num2str(GOOD_SUB(i)), '_seeds_', mm, 'mm.txt'], valid_seeds, 'precision', 20)
    dlmwrite(['./data/tractography/SUBJECTS/', num2str(GOOD_SUB(i)), '/', num2str(GOOD_SUB(i)), '_ALLGENEseeds_', mm, 'mm.txt'], allGeneSeeds, 'precision', 20)
end

% Save the list of valid subjects
dlmwrite('VALIDSEED_UnrelatedSubs.txt', GOOD_SUB, 'precision', 6)
