Mdisease=readtable('SupplementaryTable5.xlsx');
Ldisease=readtable('SupplementaryTable5.xlsx','Sheet',2);

bar_cmap = brewermap(256,'RdBu');

Minds = fliplr([7 6 9 5 1 4]);
Linds = [6 1 2];

diseases = [Mdisease.description(Minds);Ldisease.description(Linds)];

Enrich = zeros(9,2);
Enrich(1:6,1) = Mdisease.enrichmentRatio(Minds)*-1;
Enrich(7:9,2) = Ldisease.enrichmentRatio(Linds);

b = barh(Enrich,'stacked','FaceColor','flat');
 
b(1).FaceColor = bar_cmap(224,:);
b(2).FaceColor = bar_cmap(32,:);

yticks(1:length(Enrich))
yticklabels(diseases)

xlabel('Medial-lateral enrichment ratio')
set(gca,'FontSize',20)

ylim([.2 9.8])
ytickangle(0)