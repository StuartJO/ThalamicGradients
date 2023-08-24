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
not_in_AHBA_ind = ~ismember(BrainGene, AHBAgenes);

% Extract the genes from BrainGene that are not present in AHBAgenes
not_in_AHBA = BrainGene(not_in_AHBA_ind);

% Extract entrez IDs from the AHBA table
EntrezIDs = AHBA.geneEntrez_id;

% Extract entrez IDs corresponding to the indices found earlier
IDS = EntrezIDs(idx);

% Find unique entrez IDs and their corresponding indices
[UniqueIDS, UniqueIDSinds] = unique(IDS);

% Extract AHBA gene symbols corresponding to the indices found earlier
AHBAgenesBrain = AHBAgenes(idx);

%% VERY IMPORTANT

% Not all genes could be successfully mapped from their gene name in the
% Burt data to an Entrez ID. This seems to have been caused by names of
% genes being updated over time, or just I tried to match the missing ones 
% as best I could. If the gene was in http://www.meduniwien.ac.at/neuroimaging/mRNA.html,
% I used the ID it gave, otherwise I looked through https://www.ncbi.nlm.nih.gov/gene/
% to see if it had been updated
% Burt Gene      Replace gene    Entrez ID
% DMRTC1B         -- DMRTC1    -- 63947
% RP11-195F19.5   -- LOC730098 -- 730098
% PRY             -- PRY2      -- 442862
% RP11-566K11.2   -- TUBB3     -- 10381
% XAGE1B          -- XAGE1B    -- 653219

% Extract unique AHBA gene symbols corresponding to unique entrez IDs
AHBAgenesBrainUnique = AHBAgenesBrain(UniqueIDSinds);

UniqueIDS_ = [UniqueIDS;63947; 730098;442862;10381;653219];

AHBAgenesBrainUnique_ = [AHBAgenesBrainUnique;{'DMRTC1'};{'LOC730098'};{'PRY2'};{'TUBB3'};{'XAGE1B'}];

% Write the unique entrez IDs to a text file
writematrix(UniqueIDS_, './data/preprocessed/AHBAEntrez.txt')

% Write the unique AHBA gene symbols to a text file
writecell(AHBAgenesBrainUnique_, './data/preprocessed/AHBAgeneSymbol.txt')
