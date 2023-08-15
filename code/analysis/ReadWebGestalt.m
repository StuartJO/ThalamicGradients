function [DiseaseEnrich, DiseaseGeneTable, DiseaseEnrichResults] = ReadWebGestalt(FOLDER)

% This function reads and processes enrichment results from WebGestalt analysis.

% Inputs:
%   FOLDER: The folder path containing the WebGestalt enrichment result files.

% Outputs:
%   DiseaseEnrich: A table containing disease enrichment results, including columns like description, enrichmentRatio, pValue, and FDR.
%   DiseaseGeneTable: A table where each row represents a disease and columns contain the genes associated with each disease.
%   DiseaseEnrichResults: A table containing raw enrichment results including additional information.

% Find enrichment result file in the specified folder
Enrichment_file = dir([FOLDER,'/enrichment_results*']);

if ~isempty(Enrichment_file)
    
    % Read the enrichment results table from the file
    DiseaseEnrichResults = readtable([Enrichment_file.folder,'/',Enrichment_file.name],'Delimiter','tab');
    
    % Count the number of diseases included in the enrichment results
    Ndisease = size(DiseaseEnrichResults,1);
    
    % Extract genes associated with each disease
    inGenes = cell(1,Ndisease);
    for i = 1:Ndisease
        genes = [';',DiseaseEnrichResults.userId{i},';'];
        Loc = strfind(genes,';');
        for j = 1:length(Loc)-1
            inGenes{i}{j} = genes(Loc(j)+1:Loc(j+1)-1);
        end
    end
    
    % Count the number of genes for each disease
    Ngenes = zeros(Ndisease,1);
    for i = 1:Ndisease
        Ngenes(i) = length(inGenes{i});
    end
    
    % Pad the gene lists to have the same length
    MaxGenes = max(Ngenes);
    for i = 1:Ndisease
        if Ngenes(i) < MaxGenes
            inGenes{i}(Ngenes(i)+1:MaxGenes) = {''};
        end
    end
    
    % Create a table with padded gene lists and use disease descriptions as column names
    DiseaseGenesPadded = cell(MaxGenes,Ndisease);
    for i = 1:Ndisease
        DiseaseGenesPadded(:,i) = inGenes{i};
    end
    
    % Extract specific columns from DiseaseEnrichResults
    Vars2Extract={'description','enrichmentRatio','pValue','FDR'};
    colIdx = zeros(1,length(Vars2Extract));
    for i = 1:length(Vars2Extract)
        colIdx(i) = find(strcmp(DiseaseEnrichResults.Properties.VariableNames, Vars2Extract{i}));
    end
    
    % Handle duplicate disease descriptions
    DiseaseEnrichResultsDescription = DiseaseEnrichResults.description;
    [uniqueStrings, ~, stringIndex] = unique(DiseaseEnrichResultsDescription);
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
    
    % Create the DiseaseGeneTable and DiseaseEnrich tables
    DiseaseGeneTable = cell2table(DiseaseGenesPadded,'VariableNames',DiseaseEnrichResultsDescription);
    DiseaseEnrich = DiseaseEnrichResults(:,colIdx);

else
    % Return empty tables if no enrichment result file is found
    DiseaseGeneTable = table([]);
    Vars2Extract={'description','enrichmentRatio','pValue','FDR'};
    varTypes = {'string', 'double', 'double','double'};
    DiseaseEnrich = table('Size',[0 length(Vars2Extract)], 'VariableNames', Vars2Extract,'VariableTypes', varTypes);
    DiseaseEnrichResults = table([]);
end
