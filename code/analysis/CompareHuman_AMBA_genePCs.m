TractGeneNorm = load('./data/processed/TractGeneNorm.mat');

GeneNames_human = TractGeneNorm.GeneNames;

PhillipsGeneTable = readtable('PhillipsMouseThalGenes.xlsx');
GeneNames_phillips = PhillipsGeneTable.GeneSymbol;
GeneNames_phillips_upper = upper(GeneNames_phillips);

Phillips_pc1_gene = zscore(PhillipsGeneTable.PC1);

[~,idx] = ismember(GeneNames_human,upper(mouse_GeneNames));
AMBA2Human_homolog = idx(idx~=0);
Human2AMBA_homolog = find(idx~=0);
iter = 1;
annotlabels = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};

for i = 1:3
    for j = 1:3
    %subplot(3,3,sub2ind([3 3],j,i))
    figure
x = zscore(pcs_gene(Human2AMBA_homolog,i));
y = zscore(mouse_pcs_gene(AMBA2Human_homolog,j));

s = scatterfit(x,y,36,lines(1),[],2,1);
xlabel(['Human gene PC',num2str(i),' loading'])
ylabel({'Allen Mouse Brain',['gene PC',num2str(j),' loading']})
a = annotation('textbox',[0 .89 .05 .13],'String',annotlabels{iter},'FontSize',32,'EdgeColor','none');
print(['./figure_outputs/MouseVsHumanGene_',num2str(iter),'.png'],'-dpng','-r300')
iter = iter + 1;
    end   
end