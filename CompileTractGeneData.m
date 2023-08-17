function CompileTractGeneData()

GenesKept = dlmread('GenesKept.txt');

SeedGene = zeros(1811,length(GenesKept));
GeneTractDir = '/projects/kg98/stuarto/SeedReg/seeds_1.75mm_tracts_gene';
for i = 1:length(GenesKept)
SeedGene(:,i) = dlmread([GeneTractDir,'/seeds_1.75mm_tracts_',num2str(GenesKept(i)),'_gene.txt']);
disp(num2str(i))
end

SUB = dlmread('VALIDSEED_UnrelatedSubs.txt');
seed_ind = dlmread('seeds_1.75mm_fsl_seeds2use_76subs_valid_bin.txt');
Nseed = sum(seed_ind);

SeedGene_kept = SeedGene(logical(seed_ind),:);


ThalSeed = zeros(Nseed,500,length(SUB));

for j = 1:length(SUB)
    TRACTDIR = ['/projects/kg98/stuarto/SeedReg/SUBJECTS/',num2str(SUB(j)),'/tracts_921seeds'];
for i = 1:Nseed
    ThalSeed(i,:,j) = dlmread([TRACTDIR,'/thal_seed_',num2str(i),'_cp']);
end
end

ThalSeedAvg = mean(ThalSeed,3);

norm = BF_NormalizeMatrix(ThalSeedAvg(:,1:250),'scaledSigmoid');

save('Sub76_ThalData.mat','ThalSeed','ThalSeedAvg','norm','SUB','seed_ind','SeedGene_kept','GenesKept')