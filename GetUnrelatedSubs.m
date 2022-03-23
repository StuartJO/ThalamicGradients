TwinSubs = dlmread('TwinSubs.txt');
AllSubs = dlmread('HCP_SUBJECTS.txt');
IsRelated = ismember(AllSubs,TwinSubs);

UnrelatedSubs = AllSubs(~IsRelated);

[~,I] = sort(rand(1,length(UnrelatedSubs)));

UnrelatedSubs2Use = UnrelatedSubs(I(1:100));

dlmwrite('UnrelatedSubs.txt',UnrelatedSubs2Use,'Precision',6)

dlmwrite('UnrelatedAndTwinSubs.txt',[UnrelatedSubs2Use;TwinSubs],'Precision',6)