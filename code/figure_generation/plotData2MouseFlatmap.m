function plotData2MouseFlatmap(data,cmap)

cort_pc_color = MapData2Colors(data,cmap);

load('./data/ancillary/Mouseflatmap.mat','Mouseflatmap_boundary','Mouseflatmap_boundary_id')

for i = 1:38
    IND = find(Mouseflatmap_boundary_id==i);
    fill(Mouseflatmap_boundary{IND}(:,1),Mouseflatmap_boundary{IND}(:,2),cort_pc_color(i,:),'EdgeColor',[0 0 0]);
    hold on
end

NaN_regions = find(Mouseflatmap_boundary_id==0);
for i = 1:length(NaN_regions)
    IND = NaN_regions(i);
    fill(Mouseflatmap_boundary{IND}(:,1),Mouseflatmap_boundary{IND}(:,2),[.5 .5 .5],'EdgeColor',[0 0 0]);
    hold on
end

axis image
axis off