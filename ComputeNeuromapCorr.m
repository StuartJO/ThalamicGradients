addpath C:/Users/Stuart/Documents/MATLAB/gifti-1.6/

files = dir('./fsaverage164k_annots/*.gii');
for i = 1:length(files)
g= gifti(files(i).name);
neuromap(:,i) = g.cdata;
end

load('fsaverage_surface_data.mat')

for i = 1:250
neuromap_parc(i,:) = mean(neuromap(lh_rand500==i,:),1);
end


C = corr(coeff{3}(1:250,1),neuromap_parc,'Rows','complete');

% C = corr(coeff{3}(1:250,1),neuromap_parc,'Rows','complete','Type','Spearman');

for j = 1:6
    figure
    for i = 1:12
        subplot(3,4,i)
        scatter(coeff{3}(1:250,1),neuromap_parc(:,i+(12*(j-1))));
        title(files(i+(12*(j-1))).name)
    end
end

Import_receptor2tracer
tracer=receptor2tracer_map(:,2);
receptor=receptor2tracer_map(:,1);

for i = 1:length(files)
name = files(i).name;
Underloc = strfind(name,'_');
NeuroMapName = name(Underloc(1)+1:Underloc(2)-1);
Index = find(contains(tracer,NeuroMapName));
if ~isempty(Index)
NeuroMapNames{i} = receptor{Index(1)};
else
    NeuroMapNames{i} = NeuroMapName;
end
end


Chk = find(abs(C)>.4);
for i = 1:length(Chk)
    subplot(3,5,i)
    scatter(coeff{3}(1:250,1),neuromap_parc(:,Chk(i)));
end
