function MappedColors = MapData2Colors(data,cmap,limits)

if nargin < 3
    limits = [min(data) max(data)];
end

if size(data,2) == 1
    data_ = [data; limits'];
else
    data_ = [data limits];
end

 color_ind = ceil(rescale(data_,1,size(cmap,1)));

 color_ind(isnan(color_ind)) = 1;
 
 MappedColors = cmap(color_ind(1:end-2),:);