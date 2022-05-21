function vals_color = value2color(vals,cmap)

 color_ind = ceil(rescale(vals,1,size(cmap,1)));
 vals_color = cmap(color_ind,:);