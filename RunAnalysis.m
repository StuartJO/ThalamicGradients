
addpath(genpath('./'))

TractGeneData = load('CompiledTractGeneData.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);

save('./data/processed/TractGeneNorm.mat','-struct','CompiledTractData')

TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];

Ncorts = size(CompiledTractData.TractData_norm,2);

decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
decomp.GeneNames_human = CompiledTractData.GeneNames_human;
save('./data/processed/main_decomp.mat','-struct','decomp')

for i = 1:3
writematrix(decomp.score(:,i),['./data/processed/PC',num2str(i),'_thal.txt'],'Delimiter',' ')
end

load('./data/preprocessed/MNI_Seed_voxelData.mat','seed_vox_coords')

SeedDists = squareform(pdist(seed_vox_coords(logical(TractGeneData.seed_ind),:)));
writematrix(SeedDists,'./data/processed/SeedDists.txt','Delimiter',' ')

%%
TractGeneData = load('CompiledTractGeneData_AllGeneSeed.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);
save('./data/processed/TractGeneNorm_AllGeneSeed.mat','-struct','CompiledTractData')
TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];
Ncorts = size(CompiledTractData.TractData_norm,2);
decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
decomp.GeneNames_human = CompiledTractData.GeneNames_human;
save('./data/processed/decomp_AllGeneSeed.mat','-struct','decomp')

TractGeneData = load('CompiledTractGeneData_Scha400.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);
save('./data/processed/TractGeneNorm_Scha400.mat','-struct','CompiledTractData')
TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];
Ncorts = size(CompiledTractData.TractData_norm,2);
decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
decomp.GeneNames_human = CompiledTractData.GeneNames_human;
save('./data/processed/decomp_Scha400.mat','-struct','decomp')

%% Run mouse PCA

RunAMBAdecomp

%% Run alternative decompositions

[AltEmbedding,DecompName] = RunAltDecomps();

save('./data/processed/alt_decomps','AltEmbedding','DecompName')

%% Get 

main_decomp = load('main_decomp.mat');
main_genedata = load('TractGeneNorm.mat');
for i = 1:3
surrogates = csvread(['thal_surrogates_PC',num2str(i),'.csv']);  
pc_gene_coeffs = main_decomp.coeff(251:end,i);
AutoCorrReslt = PC_spatial_autocorrelation(main_decomp.score(:,i),pc_gene_coeffs,main_genedata,surrogates');

writecell([AutoCorrReslt.MostPositive num2cell(AutoCorrReslt.PC_genes_positive(1:100))],['./data/processed/PC',num2str(i),'_HumanMostPositiveSpinTested.csv'])
writecell([AutoCorrReslt.MostNegative num2cell(AutoCorrReslt.PC_genes_negative(1:100))],['./data/processed/PC',num2str(i),'_HumanMostNegativeSpinTested.csv'])
writecell([main_decomp.GeneNames_human num2cell(pc_gene_coeffs) num2cell(AutoCorrReslt.p_perm)],['./data/processed/PC',num2str(i),'_AllHumanSpinTested.csv'])

writecell(AutoCorrReslt.MostPositive,['./data/processed/PC',num2str(i),'_HumanMostPositiveSpinTested.txt'])
writecell(AutoCorrReslt.MostNegative,['./data/processed/PC',num2str(i),'_HumanMostNegativeSpinTested.txt'])

end

%%

% Uncomment the line below to remake the cortical spin test. Will take ages
%GetCorticalSpinTestPerms()

for i = 1:3
    neuromap_corrs = GetNeuromapCorrs(main_decomp.pcs_cort(:,i));
    save(['./data/processed/NeuroMapCorrs_PC',num2str(i),'.mat'],'-struct','neuromap_corrs')
end


%% Bootstrap and CV

TractGeneData = load('CompiledTractGeneData.mat');
CompiledTractData = load('./data/processed/TractGeneNorm.mat');

TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];

nVars = size(TractData_GeneData_norm,2);
nObs = size(TractData_GeneData_norm,1);

Nboot = 1000;

coeffPC1_boot = zeros(nVars,Nboot);
scorePC1_boot = zeros(nObs,Nboot);
expl_boot = zeros(Nboot,min(nVars,nObs)-1);

h = waitbar(0,'Please wait...');

if ~exist('BootstrappedInds.mat','file')
    BootInds = zeros(76,Nboot);
    BootIndsExist = 0;
else
    load('BootstrappedInds.mat','BootInds')
    BootIndsExist = 1;
end

for i = 1:Nboot   
    if ~BootIndsExist
        BootInds(:,i) = randi(76,1,76);    
    end
    TractData_norm_ = BF_NormalizeMatrix(mean(TractGeneData.ThalSeed(:,:,BootInds(:,i)),3),'scaledSigmoid');  
    TractData_GeneData_norm_ = [TractData_norm_ CompiledTractData.GeneData_norm];
    [coeff_,score_,~,~,expl_] = pca(TractData_GeneData_norm_);   
    coeffPC1_boot(:,i) = coeff_(:,1);
    scorePC1_boot(:,i) = score_(:,1); 
    expl_boot(i,:) = expl_;
    waitbar(i/nObs,h,['Finished bootstrap ',num2str(i)])
end

if ~exist('BootstrappedInds.mat','file')

save('BootstrappedInds.mat','BootInds')

end

coeffPC1_boot_tract = zeros(250,Nboot);
scorePC1_boot_tract = zeros(nObs,Nboot);
expl_boot_tract = zeros(Nboot,250);

for i = 1:Nboot 
    TractData_norm_ = BF_NormalizeMatrix(mean(TractGeneData.ThalSeed(:,1:250,BootInds(:,i)),3),'scaledSigmoid');   
    [coeff_,score_,~,~,expl_] = pca(TractData_norm_);   
    coeffPC1_boot_tract(:,i) = coeff_(:,1);
    scorePC1_boot_tract(:,i) = score_(:,1); 
    expl_boot_tract(i,:) = expl_;  
    waitbar(i/nObs,h,['Finished bootstrap ',num2str(i)])
end

tract_decomp = RunPCADecomp(CompiledTractData.TractData_norm,CompiledTractData.seed_ind,1:250,[]);
save('./data/processed/decomp_TractOnly.mat','-struct','tract_decomp')

save('./data/processed/Bootstrap_results.mat','coeffPC1_boot','scorePC1_boot','expl_boot','coeffPC1_boot_tract','scorePC1_boot_tract','expl_boot_tract')
%% Leave-one-out

h = waitbar(0,'Please wait...');

TractGeneData = load('CompiledTractGeneData.mat');

nObs = sum(TractGeneData.seed_ind);
scoresCV = nan(nObs);

SeedCVind = 1:nObs;

for i = 1:nObs
    
    SeedCVind_ = SeedCVind;
    SeedCVind_(i) = [];
    
    TractData_norm_ = BF_NormalizeMatrix(TractGeneData.ThalSeedAvg(SeedCVind_,1:250),'scaledSigmoid');
    GeneData_norm_ = BF_NormalizeMatrix(TractGeneData.SeedGene_kept(SeedCVind_,:),'scaledSigmoid');
    
    TractData_GeneData_norm_ = [TractData_norm_ GeneData_norm_];

	[coeff_,score_] = pca(TractData_GeneData_norm_);
    
    scoresCV(SeedCVind_,i) = score_(:,1);
    waitbar(i/nObs,h,['Finished CV ',num2str(i)])
end

RMSE = sqrt(nanmean((main_decomp.score(:,1) - scoresCV).^2)); % Root Mean Squared Error

save('./data/processed/loocv_result.mat','RMSE','scoresCV')