
TractGeneData = load('CompiledTractGeneData.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);

save('./data/processed/TractGeneNorm.mat','-struct','CompiledTractData')

TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];

Ncorts = size(CompiledTractData.TractData_norm,2);

decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));

save('./data/processed/main_decomp.mat','-struct','decomp')

for i = 1:3
writematrix(decomp.score(:,i),['./data/processed/PC',num2str(i),'_thal.txt'],'Delimiter',' ')
end

SeedDists = squareform(pdist(seed_vox_coords(logical(seed_ind),:)));
writematrix(SeedDists,'SeedDists.txt','Delimiter',' ')


%%
TractGeneData = load('CompiledTractGeneData_AllGeneSeed.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);
save('./data/processed/TractGeneNorm_AllGeneSeed.mat','-struct','CompiledTractData')
TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];
Ncorts = size(CompiledTractData.TractData_norm,2);
decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
save('./data/processed/decomp_AllGeneSeed.mat','-struct','decomp')

TractGeneData = load('CompiledTractGeneData_Scha400.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);
save('./data/processed/TractGeneNorm_Scha400.mat','-struct','CompiledTractData')
TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];
Ncorts = size(CompiledTractData.TractData_norm,2);
decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
save('./data/processed/decomp_Scha400.mat','-struct','decomp')

%% Run mouse PCA

RunAMBAdecomp

%% Run alternative decompositions

[AltEmbedding,DecompName] = RunAltDecomps();

save('./data/processed/alt_decomps','AltEmbedding','DecompName')

%% Get 

main_decomp = load('main_decomp.mat');

for i = 1:3
surrogates = csvread(['thal_surrogates_PC',num2str(i),'.csv']);  
pc_gene_coeffs = main_decomp.coeff(251:end,i);
AutoCorrReslt = PC_spatial_autocorrelation(main_decomp.score(:,i),pc_gene_coeffs,surrogates');

writecell([AutoCorrReslt.MostPositive num2cell(AutoCorrReslt.PC_genes_positive(1:100))],['./data/processed/PC',num2str(i),'_HumanMostPositiveSpinTested.csv'])
writecell([AutoCorrReslt.MostNegative num2cell(AutoCorrReslt.PC_genes_negative(1:100))],['./data/processed/PC',num2str(i),'_HumanMostNegativeSpinTested.csv'])
writecell([main_decomp.GeneNames_human num2cell(pc_gene_coeffs) num2cell(AutoCorrReslt.p_perm)],['./data/processed/PC',num2str(i),'_AllHumanSpinTested.csv'])

writecell(AutoCorrReslt.MostPositive,['./data/processed/PC',num2str(i),'_HumanMostPositiveSpinTested.txt'])
writecell(AutoCorrReslt.MostNegative,['./data/processed/PC',num2str(i),'_HumanMostNegativeSpinTested.txt'])

end

%%


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

main_decomp = load('main_decomp.mat');

figure('Position',[146 301 1073 590])
offset_x = -.01;
offset_y = .02;

subplot(2,3,2)
histogram(corr(main_decomp.score(:,1),scorePC1_boot))
xlabel('Bootstrapped score correlation with original')
ylabel('Count')
title('Gene+tract')
addPlotLabel('b',gca,24,[offset_x offset_y])
subplot(2,3,3)
histogram(corr(main_decomp.coeff(:,1),coeffPC1_boot))
xlabel('Bootstrapped loading correlation with original')
ylabel('Count')
title('Gene+tract')
addPlotLabel('c',gca,24,[offset_x offset_y])
subplot(2,3,5)
histogram(corr(tract_decomp.score(:,1),scorePC1_boot_tract))
xlabel('Bootstrapped score correlation with original')
ylabel('Count')
title('Tract')
addPlotLabel('e',gca,24,[offset_x offset_y])
subplot(2,3,6)
histogram(corr(tract_decomp.coeff(:,1),coeffPC1_boot_tract))
xlabel('Bootstrapped loading correlation with original')
ylabel('Count')
title('Tract')
addPlotLabel('f',gca,24,[offset_x offset_y])

subplot(2,3,1)
for i = 1:5
V(i) = Violin(expl_boot(:,i), i,'ShowMean',true,'ShowData',false,'ViolinAlpha',.25);
V(i).ViolinColor = [.5 .5 .5];
V(i).MeanPlot.Color = [0 0 0];
V(i).EdgeColor = [0 0 0];
V(i).ViolinPlot.LineWidth = .5;
V(i).BoxColor = [0 0 0];
end
hold on
plot(main_decomp.explained(1:5),'r','LineWidth',2)
xlabel('PC')
ylabel('Variance explained')
title('Gene+tract')
addPlotLabel('a',gca,24,[offset_x offset_y])
subplot(2,3,4)
for i = 1:5
V(i) = Violin(expl_boot_tract(:,i), i,'ShowMean',true,'ShowData',false,'ViolinAlpha',.25);
V(i).ViolinColor = [.5 .5 .5];
V(i).MeanPlot.Color = [0 0 0];
V(i).EdgeColor = [0 0 0];
V(i).ViolinPlot.LineWidth = .5;
V(i).BoxColor = [0 0 0];
end
hold on
plot(tract_decomp.explained(1:5),'r','LineWidth',2)
xlabel('PC')
title('Tract')
ylabel('Variance explained')
addPlotLabel('d',gca,24,[offset_x offset_y])

exportgraphics(gcf,'bootstrap_result.png','Resolution',300)

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

%histogram(nanvar(scoresCV,[],2))

RMSE = sqrt(nanmean((main_decomp.score(:,1) - scoresCV).^2)); % Root Mean Squared Error

figure
histogram(RMSE)

xlabel('RMSE of PC1 scores')
ylabel({'Number of','leave-one-out iterations'})
set(gca,'FontSize',20)

exportgraphics(gcf,'loocv_result.png','Resolution',300)