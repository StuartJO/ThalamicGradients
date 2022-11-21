NeuronSubtypes={'Habenula_Tac2','Rora','Gad2-Ahi1'};

for i = 1:3
    intable=readtable(['NegativeNeuron_',NeuronSubtypes{i},'Genes.csv']);
    inGenes{i}=intable.Gene;
end
for i = 1:3
    intable=readtable(['PositiveNeuron_',NeuronSubtypes{i},'Genes.csv']);
    inGenes{i+3}=intable.Gene;
end

for i = 1:6
    Ngenes(i) = length(inGenes{i});
end

MaxGenes = max(Ngenes);

for i = 1:6
    if Ngenes(i) < MaxGenes
    inGenes{i}(Ngenes(i)+1:MaxGenes) = {''};
    end
end

for i = 1:6
    SubtybeGenesPadded(:,i) = inGenes{i};
end

for i = 1:3
NeuronSubtypeGeneVarName{i} = ['Medial_',NeuronSubtypes{i}];
NeuronSubtypeGeneVarName{i+3} = ['Lateral_',NeuronSubtypes{i}];
end

NeuronSubtypeGenes = cell2table(SubtybeGenesPadded,'VariableNames',NeuronSubtypeGeneVarName);

writetable(NeuronSubtypeGenes,'NeuronSubtypeGenes.xlsx')