%% RunAnalysis_2

addpath(genpath('./'))

% Uncomment the line below to remake the cortical spin test. Will take ages
% GetCorticalSpinTestPerms()

decomp_rand500 = load('./data/processed/decomp_rand500.mat');

for i = 1:3
    neuromap_corrs = GetNeuromapCorrs(decomp_rand500.pcs_cort(:,i));
    save(['./data/processed/NeuroMapCorrs_PC',num2str(i),'.mat'],'-struct','neuromap_corrs')
end


%% Bootstrap and CV

TractGeneData = load('./data/preprocessed/CompiledTractGeneData_rand500.mat');
CompiledTractData = load('./data/processed/TractGeneNorm_rand500.mat');

TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];

nVars = size(TractData_GeneData_norm,2);
nObs = size(TractData_GeneData_norm,1);

Nboot = 1000;

coeffPC1_boot = zeros(nVars,Nboot);
scorePC1_boot = zeros(nObs,Nboot);
expl_boot = zeros(Nboot,min(nVars,nObs)-1);

h = waitbar(0,'Please wait...');

if ~exist('./data/processed/BootstrappedInds.mat','file')
    BootInds = zeros(76,Nboot);
    BootIndsExist = 0;
else
    load('./data/processed/BootstrappedInds.mat','BootInds')
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

if ~exist('./data/processed/BootstrappedInds.mat','file')

save('./data/processed/BootstrappedInds.mat','BootInds')

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

TractGeneData = load('./data/preprocessed/CompiledTractGeneData_rand500.mat');

nObs = sum(TractGeneData.seed_ind);
scoresCV = nan(nObs);

SeedCVind = 1:nObs;

for i = 1:nObs
    
    SeedCVind_ = SeedCVind;
    SeedCVind_(i) = [];
    
    TractData_norm_ = BF_NormalizeMatrix(TractGeneData.ThalSeedAvg(SeedCVind_,1:250),'scaledSigmoid');
    GeneData_norm_ = BF_NormalizeMatrix(TractGeneData.ThalSeedGenesKept(SeedCVind_,:),'scaledSigmoid');
    
    TractData_GeneData_norm_ = [TractData_norm_ GeneData_norm_];

	[coeff_,score_] = pca(TractData_GeneData_norm_);
    
    scoresCV(SeedCVind_,i) = score_(:,1);
    waitbar(i/nObs,h,['Finished CV ',num2str(i)])
end

RMSE = sqrt(nanmean((decomp_rand500.score(:,1) - scoresCV).^2)); % Root Mean Squared Error

save('./data/processed/loocv_result.mat','RMSE','scoresCV')