function x = shellsort(x)
%--------------------------------------------------------------------------
% Syntax:       sx = shellsort(x);
%               
% Inputs:       x is a vector of length n
%               
% Outputs:      sx is the sorted (ascending) version of x
%               
% Description:  This function sorts the input array x in ascending order
%               using the Shell sort algorithm
%               
% Complexity:   *Depends on gap sequence*    best-case performance
%               *Depends on gap sequence*    average-case performance
%               *Depends on gap sequence*    worst-case performance
%               O(1)                         auxiliary space
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 5, 2014
%--------------------------------------------------------------------------
% Compute gap sequence
% Note: In practice, one should preallocate memory for gaps
n = length(x);
gap = 1;
gaps = [];
while (ceil(gap) <= n)
    gaps = [gaps ceil(gap)]; %#ok
    gap = 2.25 * gap + 1; % Tokuda's gap sequence
end
Ngaps = length(gaps);
% Insertion sort
for k = Ngaps:-1:1
    % Perform insertion sort with gap = gaps(k)
    x = insertionsorti(x,gaps(k),n);    
end
end
function x = insertionsorti(x,gap,n)
% Perform insertion with given gap
% Note: In practice, x xhould be passed by reference
% Performs insertion sort with given gap
for j = (gap + 1):n
    pivot = x(j);
    i = j;
    while ((i > gap) && (x(i - gap) > pivot))
        x(i) = x(i - gap);
        i = i - gap;
    end
    x(i) = pivot;
    imagesc(x)
    pause(.1)
end
end
