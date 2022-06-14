
pc_thal = (zscore(score{3}(:,1)'));
[pc_thal_sort, pc_thal_sort_ind] = sort(zscore(score{3}(:,1)'),'descend');

imagesc(pc_thal)
colormap(turbo(256))
array = pc_thal;
                imagesc(array)
                axis off
            pause(.1)
print(['./GIF/SORT_PCTHAL_',num2str(1),'.png'],'-dpng')
ITER = 2;
for i = 1:length(array)
    % Pass up to the last un-sorted element.
    for j = 1:length(array)-i
        % If elements are in the wrong order, swap them.
        if(array(j)<array(j+1))
            temp = array(j+1);
            array(j+1) = array(j);
            array(j) = temp;

        end
    end
    array_diff(i) = sum((array-pc_thal_sort).^2);
    if array_diff(i) == 0
        imagesc(array)
        axis off
        pause(.1)
        print(['./GIF/SORT_PCTHAL_',num2str(ITER),'.png'],'-dpng')  
        break
    end
    
    if mod(i,3) == 0
                imagesc(array)
                axis off
            pause(.1)
print(['./GIF/SORT_PCTHAL_',num2str(ITER),'.png'],'-dpng')    
ITER = ITER + 1;
    end
end