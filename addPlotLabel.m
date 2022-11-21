function addPlotLabel(TEXT,AXIS,FONTSIZE,OFFSET)

if nargin < 2
    AXIS = gca;
end

if nargin < 3
    FONTSIZE = 24;
end

if nargin < 4
    OFFSET = [0 0];
end

axis_pos = get(AXIS,'OuterPosition');
annotYpos = find_point_on_line(axis_pos(2),axis_pos(2)+axis_pos(4),.9);
annotation('textbox',[axis_pos(1)+OFFSET(1) annotYpos+OFFSET(2) axis_pos(3)*.1 axis_pos(2)+axis_pos(4)-annotYpos],'String',TEXT,'FontSize',FONTSIZE,'EdgeColor','none')
