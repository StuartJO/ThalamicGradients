Medial_folder='Project_wg_result1665153140';
Lateral_folder='Project_wg_result1665153115';


for MedLat = 1:2

    if MedLat == 1
        Enrichment_files = dir([Medial_folder,'/enrichment_results*']);
    else
        Enrichment_files = dir([Lateral_folder,'/enrichment_results*']);
    end

DiseaseEnrichResults = importDiseaseEnrichfile([Enrichment_files.folder,'/',Enrichment_files.name]);

Ndisease = size(DiseaseEnrichResults,1);
inGenes = cell(1,Ndisease);
for i = 1:Ndisease
    genes = [';',DiseaseEnrichResults.enriched_genes{i},';'];
    Loc = strfind(genes,';');
    for j = 1:length(Loc)-1
       inGenes{i}{j} = genes(Loc(j)+1:Loc(j+1)-1); 
    end
end

Ngenes = zeros(Ndisease,1);

for i = 1:Ndisease
    Ngenes(i) = length(inGenes{i});
end

MaxGenes = max(Ngenes);

for i = 1:Ndisease
    if Ngenes(i) < MaxGenes
    inGenes{i}(Ngenes(i)+1:MaxGenes) = {''};
    end
end

DiseaseGenesPadded = cell(MaxGenes,Ndisease);

for i = 1:Ndisease
    DiseaseGenesPadded(:,i) = inGenes{i};
end


if MedLat == 1
    MedDiseaseGeneTable = cell2table(DiseaseGenesPadded,'VariableNames',DiseaseEnrichResults.description);
    MedDiseaseEnrich = DiseaseEnrichResults(:,[2 7 8 9]);
else
    LatDiseaseGeneTable = cell2table(DiseaseGenesPadded,'VariableNames',DiseaseEnrichResults.description);
    LatDiseaseEnrich = DiseaseEnrichResults(:,[2 7 8 9]);
end

end

writetable(MedDiseaseEnrich,'DiseaseEnrichment.xlsx','WriteMode','overwrite','Sheet','MedialDiseaseEnrich')

writetable(LatDiseaseEnrich,'DiseaseEnrichment.xlsx','WriteMode','Append','Sheet','LateralDiseaseEnrich')

writetable(MedDiseaseGeneTable,'DiseaseEnrichment.xlsx','WriteMode','Append','Sheet','MedialDiseaseGenes')

writetable(LatDiseaseGeneTable,'DiseaseEnrichment.xlsx','WriteMode','Append','Sheet','LateralDiseaseGenes')