function plot_gene_trajectories(PC,savename)

if nargin < 1
    PC = 1;
end

if nargin < 2
    savename = [];
end

bar_cmap = brewermap(256,'RdBu');
size_cmap = size(bar_cmap,1);
 pos_bar_cmap = flipud(bar_cmap(1:128,:));
 neg_bar_cmap = (bar_cmap(129:size_cmap,:));
RedColor = bar_cmap(32,:);
BlueColor = bar_cmap(224,:);
figure('Position',[68 192 1317 641])
for i = 1:2
subplot(1,2,i)
if i == 1
    
   traj = readtable(['negative_trajectories_PC',num2str(PC),'.csv']); 
   plotcmap = neg_bar_cmap;  
   Plotlabel = 'a';
else
    
   traj = readtable(['positive_trajectories_PC',num2str(PC),'.csv']); 
   plotcmap = pos_bar_cmap;
   Plotlabel = 'b';
end


GenestrajID = unique(traj.symbol);

enrichlimits = [min(traj.fit) max(traj.fit)];

hold on
for i = 1:length(GenestrajID)
   IND = strcmp(traj.symbol,GenestrajID{i});
   x = traj.age(IND);
   y = traj.fit(IND);
   
   plotcolor = MapData2Colors(max(y),plotcmap,enrichlimits);
   
   pos_traj_plot(i) = plot(x,y,'Color',[plotcolor 1]);
end

set(gca,'XScale','log')
ylimits = ylim;
ylim([ylimits(1) find_point_on_line(ylimits(1),ylimits(2),1.1)])
ylimits = ylim;
plot([280 280],ylimits,'Color','k')
t1 = text(280,find_point_on_line(ylimits(1),ylimits(2),.95),'  postnatal','HorizontalAlignment','left','FontSize',20);
t2 = text(280,find_point_on_line(ylimits(1),ylimits(2),.95),'prenatal  ','HorizontalAlignment','right','FontSize',20);
xlabel('Age in days (log)')
ylabel('Relative gene expression')
title(['PC',num2str(PC)])
set(gca,'FontSize',20)
axis_pos = get(gca,'OuterPosition');
    annotYpos = find_point_on_line(axis_pos(2),axis_pos(2)+axis_pos(4),.9);
    annotation('textbox',[axis_pos(1)+.015 annotYpos axis_pos(3).*.1 axis_pos(2)+axis_pos(4)-annotYpos],'String',Plotlabel,'FontSize',24,'EdgeColor','none')
end

if ~isempty(savename)
    exportgraphics(gcf,savename,'Resolution',300)
end