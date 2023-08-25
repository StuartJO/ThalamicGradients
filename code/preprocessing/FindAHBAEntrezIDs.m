function FindAHBAEntrezIDs()

% Read the table containing AHBA (Allen Human Brain Atlas) entrez IDs from an Excel file
AHBA = readtable('./data/preprocessed/AHBA_entrez_ids.xlsx');

% Read the table containing brain genes information from an Excel file
BrainGeneTable = readtable('./data/preprocessed/BrainGenes.xlsx');

% Extract the 'BrainGenes' column from the BrainGeneTable, which likely contains gene symbols
BrainGene = BrainGeneTable.BrainGenes;

% Extract gene symbols from the AHBA table
AHBAgenes = AHBA.geneSymbol;

% Find the indices of AHBA genes that are present in the BrainGene list
idx = find(ismember(AHBAgenes, BrainGene));

% Find the indices of genes in BrainGene that are not present in AHBAgenes
%not_in_AHBA_ind = ~ismember(BrainGene, AHBAgenes);

% Extract the genes from BrainGene that are not present in AHBAgenes
%not_in_AHBA = BrainGene(not_in_AHBA_ind);

% Extract entrez IDs from the AHBA table
EntrezIDs = AHBA.geneEntrez_id;

% Extract entrez IDs corresponding to the indices found earlier
IDS = EntrezIDs(idx);

% Find unique entrez IDs and their corresponding indices
[UniqueIDS, UniqueIDSinds] = unique(IDS);

% Extract AHBA gene symbols corresponding to the indices found earlier
AHBAgenesBrain = AHBAgenes(idx);

% Extract unique AHBA gene symbols corresponding to unique entrez IDs
AHBAgenesBrainUnique = AHBAgenesBrain(UniqueIDSinds);

% Write the unique entrez IDs to a text file
writematrix(UniqueIDS, './data/preprocessed/AHBAEntrez.txt')

% Write the unique AHBA gene symbols to a text file
writecell(AHBAgenesBrainUnique, './data/preprocessed/AHBAgeneSymbol.txt')
