function CompileTractGeneData()

GenesKept = dlmread('GenesKept.txt');

SeedGene = zeros(1811,length(GenesKept));
GeneTractDir = './data/AHBA_wholebrain/seeds_1.75mm_tracts_gene';
for i = 1:length(GenesKept)
SeedGene(:,i) = dlmread([GeneTractDir,'/seeds_1.75mm_tracts_',num2str(GenesKept(i)),'_gene.txt']);
disp(num2str(i))
end

SUB = dlmread('VALIDSEED_UnrelatedSubs.txt');
seed_ind = dlmread('./data/preprocessing/seeds_1.75mm_fsl_seeds2use_76subs_valid_bin.txt');
Nseed = sum(seed_ind);

ThalSeedGenesKept = SeedGene(logical(seed_ind),:);

ThalSeed = zeros(Nseed,500,length(SUB));

for j = 1:length(SUB)
    TRACTDIR = ['./data/tractography/SUBJECTS/',num2str(SUB(j)),'/tracts_921seeds'];
    
    for i = 1:Nseed
        ThalSeed(i,:,j) = dlmread([TRACTDIR,'/random_500/thal_seed_',num2str(i),'_cp']);
    end
    
end

ThalSeedAvgFull = mean(ThalSeed,3);

ThalSeedAvg = ThalSeedAvgFull(:,1:250);

save('./data/preprocessed/CompiledTractGeneData_Rand500.mat','ThalSeed','ThalSeedAvg','SUB','seed_ind','ThalSeedGenesKept','GenesKept')

%%

ThalSeed = zeros(Nseed,400,length(SUB));

for j = 1:length(SUB)
    TRACTDIR = ['./data/tractography/SUBJECTS/',num2str(SUB(j)),'/tracts_921seeds'];
    
    for i = 1:Nseed
        ThalSeed(i,:,j) = dlmread([TRACTDIR,'/Schaefer400_17net/thal_seed_',num2str(i),'_cp']);
    end
    
end

ThalSeedAvgFull = mean(ThalSeed,3);

ThalSeedAvg = ThalSeedAvgFull(:,1:200);

save('./data/preprocessed/CompiledTractGeneData_Scha_400.mat','ThalSeed','ThalSeedAvg','SUB','seed_ind','ThalSeedGenesKept','GenesKept')

%%

seed_ind = dlmread('./data/preprocessing/gene_masked_seeds_1.75mm_ind.txt');
Nseed = sum(seed_ind);

ThalSeedGenesKept = SeedGene(logical(seed_ind),:);

ThalSeed = zeros(Nseed,500,length(SUB));

for j = 1:length(SUB)
    TRACTDIR = ['./data/tractography/SUBJECTS/',num2str(SUB(j)),'/tracts_ALLGENEseeds'];
    
    for i = 1:Nseed
        ThalSeed(i,:,j) = dlmread([TRACTDIR,'/random_500/thal_seed_',num2str(i),'_cp']);
    end
    
end

ThalSeedAvgFull = mean(ThalSeed,3);

ThalSeedAvg = ThalSeedAvgFull(:,1:250);

save('./data/preprocessed/CompiledTractGeneData_AllGeneSeed.mat','ThalSeed','ThalSeedAvg','SUB','seed_ind','ThalSeedGenesKept','GenesKept')

