MainGeneSeed = load('./data/processed/TractGeneNorm_rand500.mat');
AllGeneSeed = load('./data/processed/TractGeneNorm_AllGeneSeed.mat');

load('./data/ancillary/MNI_Seed_voxelData.mat','seed_vox_coords','seed_mni_coords')

seedIn = MainGeneSeed.seed_ind + AllGeneSeed.seed_ind;

CoverageCmap = [1 1 1; lines(2)];

c = PlotThalGradientSlicesAlt(seedIn,seed_vox_coords,CoverageCmap,'Coverage',2.1);

caxis([-.5 2.5])

c.Ticks = 0:2;
c.TickLabels = {'Full','Gene','Used'};

print('./figure_outputs/ThalamicCoverage','-dpng','-r300')

print('./figure_outputs/FigureS5.svg','-dsvg','-r300')
