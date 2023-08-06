function plotData2MouseFlatmap(data,cmap,offset)

if nargin < 3
   offset = [0 0 1]; 
end

cort_pc_color = MapData2Colors(data,cmap);

load('./data/ancillary/Mouseflatmap.mat','Mouseflatmap_boundary','Mouseflatmap_boundary_id')

Mouseflatmap_boundary_offset = cell(42,1);

Mouseflatmap_boundary_centre = mean(Mouseflatmap_boundary{1});

for i = 1:42
    Mouseflatmap_boundary_offset{i} = ((Mouseflatmap_boundary{i}-Mouseflatmap_boundary_centre).*offset(3))+[offset(1:2) 0];
end

for i = 1:38
    IND = find(Mouseflatmap_boundary_id==i);
    fill(Mouseflatmap_boundary_offset{IND}(:,1),Mouseflatmap_boundary_offset{IND}(:,2),cort_pc_color(i,:),'EdgeColor',[0 0 0]);
    hold on
end

NaN_regions = find(Mouseflatmap_boundary_id==0);
for i = 1:length(NaN_regions)
    IND = NaN_regions(i);
    fill(Mouseflatmap_boundary_offset{IND}(:,1),Mouseflatmap_boundary_offset{IND}(:,2),[.5 .5 .5],'EdgeColor',[0 0 0]);
    hold on
end

%axis image
%axis off