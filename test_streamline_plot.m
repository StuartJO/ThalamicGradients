
xlimits_side = [22.8127  156.9987];
ylimits_side = [24.5933  195.4591];
zlimits_side = [27.7807  149.4851];

t = 60;
%views = [linspace(0,90,t)' linspace(90,0,t)'];
% xlimits_range = [linspace(slice_xlimits(1),xlimits_side(1),t)' linspace(slice_xlimits(2),xlimits_side(2),t)'];
% ylimits_range = [linspace(slice_ylimits(1),ylimits_side(1),t)' linspace(slice_ylimits(2),ylimits_side(2),t)'];
% zlimits_range = [linspace(slice_zlimits(1),zlimits_side(1),t)' linspace(slice_zlimits(2),zlimits_side(2),t)'];

views = [repmat([0 90],10,1);linspace(0,90,t-10)' linspace(90,0,t-10)'];
[xlimits_range,ylimits_range,zlimits_range] = getLimitsRange([slice_xlimits; slice_ylimits; slice_zlimits]...
    ,[xlimits_side; ylimits_side; zlimits_side],t,10,[]);

Cortalphavals_rev = [0 0 0 0 0 0 0 0 0 0 linspace(0,0.25,50)];

slice_alphavals_rev = [1 1 1 1 1 1 1 1 1 1 linspace(1,0,20) zeros(1,30)];

tracks = read_mrtrix_tracks('192035inMNI_thal_seed_620_cortical.tck');
clear tracts_MNI
for i = 1:length(tracks.data)
tracts_MNI{i} = MNI_mm2vox(tracks.data{i},'mm');
end

for j = 1:t
    if exist('streamline_handle','var')
    delete(streamline_handle)
    end
    for i = 1:length(tracts_MNI)
        s = tracts_MNI{i};
        streamline_length = size(s,1);

        Sr = round(find_point_on_line(1,streamline_length,j/60));
        
        S = s(1:Sr,:)+[.5 .5 0];
        
        if size(S,1) > 1

        [~,SDIR] = max(sum(abs(diff(S))));    

        if SDIR == 1
            SCOL = [1 0 0 .2];
            col = 'r';
        elseif SDIR == 2
            SCOL = [0 1 0 .2];
            col = 'g';
        else
            SCOL = 	[0 0 1 .2];
            col = 'b';
        end


%         x = S(:,1);
%         y = S(:,2); 
%         z = S(:,3); 

        %handles.handle_plotCD(i) = plot3(x,y,z,'LineWidth',1,'Color',SCOL,'Clipping','off');
        %handles.handle_plotCD(i) = plot3t(x,y,z,.5,col);
        StepRate = 1;

        streamline_handle(i) = patch([S(:,1)' NaN],[S(:,2)' NaN],[S(:,3)' NaN],0); 
        joint = ([S(2:end,:); NaN NaN NaN]-S);
        joint(end,:) = joint(end-1,:);
        temp_joint = joint;
        joint(:,1) = temp_joint(:,1);
        joint(:,2) = temp_joint(:,2);
        cdata = [abs(joint./StepRate); NaN NaN NaN];
        cdata = reshape(cdata,length(cdata),1,3);
        set(streamline_handle(i),'CData', cdata, 'EdgeColor','interp','FaceColor','interp','Clipping','off') 

        hold on

        end
    end

slicesurf.FaceAlpha=slice_alphavals_rev(j);
pMNIsurface.FaceAlpha = Cortalphavals_rev(j);
pMNIsurfaceR.FaceAlpha = Cortalphavals_rev(j);
    view(views(j,:))
    xlim(xlimits_range(j,:));
    ylim(ylimits_range(j,:));
    zlim(zlimits_range(j,:));
    pause(.1)
end


% 
% tracks = read_mrtrix_tracks('192035inMNI_thal_seed_620_cortical.tck');
% clear tracts_MNI
% for i = 1:length(tracks.data)
% tracts_MNI{i} = MNI_mm2vox(tracks.data{i},'mm');
% end
% 
% for j = 60
%     if exist('streamline_handle','var')
%     delete(streamline_handle)
%     end
%     for i = 1:length(tracts_MNI)
%         s = tracts_MNI{i};
%         streamline_length = size(s,1);
% 
%         Sr = round(find_point_on_line(1,streamline_length,j/60));
%         
%         S = s(1:Sr,:);
%         
%         if size(S,1) > 1
% 
%         [~,SDIR] = max(sum(abs(diff(S))));    
% 
%         if SDIR == 1
%             SCOL = [1 0 0 .2];
%             col = 'r';
%         elseif SDIR == 2
%             SCOL = [0 1 0 .2];
%             col = 'g';
%         else
%             SCOL = 	[0 0 1 .2];
%             col = 'b';
%         end
% 
% 
% %         x = S(:,1);
% %         y = S(:,2); 
% %         z = S(:,3); 
% 
%         %handles.handle_plotCD(i) = plot3(x,y,z,'LineWidth',1,'Color',SCOL,'Clipping','off');
%         %handles.handle_plotCD(i) = plot3t(x,y,z,.5,col);
%         StepRate = 1;
% 
%         streamline_handle(i) = patch([S(:,1)' NaN],[S(:,2)' NaN],[S(:,3)' NaN],0); 
%         joint = ([S(2:end,:); NaN NaN NaN]-S);
%         joint(end,:) = joint(end-1,:);
%         temp_joint = joint;
%         joint(:,1) = temp_joint(:,1);
%         joint(:,2) = temp_joint(:,2);
%         cdata = [abs(joint./StepRate); NaN NaN NaN];
%         cdata = reshape(cdata,length(cdata),1,3);
%         set(streamline_handle(i),'CData', cdata, 'EdgeColor','interp','FaceColor','interp','Clipping','off') 
% 
%         hold on
% 
%         end
%     end
% pause(.1)
% end
