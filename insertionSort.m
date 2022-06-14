function A = insertionSort(A)
%% Insertion Sort Algorithm
% Insertion sort is a simple sorting algorithm that builds the final sorted
% array (or list) one item at a time. It is much less efficient on large lists 
% than more advanced algorithms such as quicksort, heapsort, or merge sort.
%
% Example: 
%
% consider we want to sort elements in the vector A
% >> A = rand(1,10);
% >> SortedA = insertionSort(A);
%
% -------------------------------------------------
% code by: Reza Ahmadzadeh (reza.ahmadzadeh@iit.it
% -------------------------------------------------
% Reference: https://en.wikipedia.org/wiki/Insertion_sort
%
%
for ii=2:length(A)
    jj = ii;
    while (jj > 1 && A(jj-1) > A(jj))
        [A(jj),A(jj-1)] = deal(A(jj-1),A(jj));
        jj = jj - 1;
    end
end
end