function SeedQC()

type = 'fsl';
mm = '1.75';

SeedThr = .7;
SubThr = .85;

SUB = dlmread('UnrelatedSubs.txt');

Nseeds = length(dlmread(['./data/preprocessed/seeds_',mm,'mm_ind.txt']));
DATA = zeros(Nseeds,10);

for i = 1:length(SUB)
    DATAtemp = dlmread(['./data/tractography/SUBJECTS/',num2str(SUB(i)),'/masked_seed_ind_',mm,'mm_',type,'_rows_',num2str(SUB(i)),'.txt']);
    DATA(DATAtemp,i) = 1;
end

InGene = dlmread('gene_masked_seeds_1.75mm_ind.txt');
InGeneLogical = false(Nseeds,1);
InGeneLogical(InGene) = true;
DATA(~InGeneLogical,:) = 0;

% imagesc(DATA)
% ylabel('Seed number')
% xlabel('Subject')
% set(gca,'FontSize',20)
% title('Yellow = in mask, blue = outside mask')

C = corr(DATA);

% figure
% imagesc(C)
% cc = colorbar;
% cc.Label.String = 'Correlation';
% set(gca,'FontSize',20)

mean_corr = mean(C);

DATA_ORIG = DATA;
DATA_KEEP_IND = mean_corr>=SeedThr;
GOOD_SUB = SUB(DATA_KEEP_IND);

DATA = DATA_ORIG(:,DATA_KEEP_IND);

SubValid = sum(DATA,1)./Nseeds;

NSubValid = length(SubValid);

SeedValid = sum(DATA,2)./length(GOOD_SUB);

ValidSeed = SeedValid>SubThr;

dlmwrite(['./data/preprocessing/seeds_',mm,'mm_',type,'_seeds2use_',num2str(NSubValid),'subs_valid_bin.txt'],double(ValidSeed)')
dlmwrite(['./data/preprocessing/seeds_',mm,'mm_',type,'_seeds2use_',num2str(NSubValid),'subs_valid_thr.txt'],SeedValid')

for i = 1:length(GOOD_SUB)

    tracks = read_mrtrix_tracks(['./data/tractography/SUBJECTS/',num2str(GOOD_SUB(i)),'/',num2str(GOOD_SUB(i)),'_seeds_',mm,'mm_fsl.tck']);
    tracts = tracks.data;
    valid_seeds = tracts(ValidSeed)';
    allGeneSeeds = tracts(InGeneLogical)';
    dlmwrite(['./data/tractography/SUBJECTS/',num2str(GOOD_SUB(i)),'/',num2str(GOOD_SUB(i)),'_seeds_',mm,'mm.txt'],valid_seeds,'precision',20)

    dlmwrite(['./data/tractography/SUBJECTS/',num2str(GOOD_SUB(i)),'/',num2str(GOOD_SUB(i)),'_ALLGENEseeds_',mm,'mm.txt'],allGeneSeeds,'precision',20)
end

dlmwrite('VALIDSEED_UnrelatedSubs.txt',GOOD_SUB,'precision',6)

