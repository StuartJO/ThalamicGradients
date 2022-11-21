


neuromap_corrs = GetNeuromapCorrs();

save('./data/processed/NeuroMapCorrs.mat','neuromap_corrs')

annotlabels = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};

NeuromapLgd = '';
for i = 1:length(SigMapDescrips)
    NeuromapLgd = [NeuromapLgd,annotlabels{i},', ',SigMapDescrips{i},'. '];
end
