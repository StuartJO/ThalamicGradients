function [DiseaseEnrich,DiseaseGeneTable,DiseaseEnrichResults] = ReadWebGestalt(FOLDER)

Enrichment_file = dir([FOLDER,'/enrichment_results*']);

if ~isempty(Enrichment_file)

DiseaseEnrichResults = readtable([Enrichment_file.folder,'/',Enrichment_file.name],'Delimiter','tab');

Ndisease = size(DiseaseEnrichResults,1);
inGenes = cell(1,Ndisease);
for i = 1:Ndisease
    genes = [';',DiseaseEnrichResults.userId{i},';'];
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

Vars2Extract={'description','enrichmentRatio','pValue','FDR'};
colIdx = zeros(1,length(Vars2Extract));
for i = 1:length(Vars2Extract)
    colIdx(i) = find(strcmp(DiseaseEnrichResults.Properties.VariableNames, Vars2Extract{i}));
end

DiseaseEnrichResultsDescription = DiseaseEnrichResults.description;

[uniqueStrings, ~, stringIndex] = unique(DiseaseEnrichResultsDescription);

% Count the number of occurrences of each unique string
stringCounts = histc(stringIndex, 1:numel(uniqueStrings));

Duplicates = find(stringCounts~=1);

if ~isempty(Duplicates)
    for i = 1:length(Duplicates)
    IN = find(ismember(DiseaseEnrichResultsDescription,uniqueStrings{Duplicates(i)}));

        for j = 2:length(IN)
            DiseaseEnrichResultsDescription{IN(j)} = [DiseaseEnrichResultsDescription{IN(j)},'_',num2str(j-1)];
        end

    end
end

DiseaseGeneTable = cell2table(DiseaseGenesPadded,'VariableNames',DiseaseEnrichResultsDescription);
DiseaseEnrich = DiseaseEnrichResults(:,colIdx);

else
    DiseaseGeneTable = table([]);
    Vars2Extract={'description','enrichmentRatio','pValue','FDR'};
    varTypes = {'string', 'double', 'double','double'};
    DiseaseEnrich = table('Size',[0 length(Vars2Extract)], 'VariableNames', Vars2Extract,'VariableTypes', varTypes);
    DiseaseEnrichResults = table([]);
end

% writetable(MedDiseaseEnrich,'DiseaseEnrichment.xlsx','WriteMode','overwrite','Sheet','MedialDiseaseEnrich')
% writetable(LatDiseaseEnrich,'DiseaseEnrichment.xlsx','WriteMode','Append','Sheet','LateralDiseaseEnrich')
% writetable(MedDiseaseGeneTable,'DiseaseEnrichment.xlsx','WriteMode','Append','Sheet','MedialDiseaseGenes')
% writetable(LatDiseaseGeneTable,'DiseaseEnrichment.xlsx','WriteMode','Append','Sheet','LateralDiseaseGenes')