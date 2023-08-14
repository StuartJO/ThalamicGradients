xlsxFiles = dir(fullfile('./SourceDataTables', '*.xlsx'));

Nfiles = length(xlsxFiles);

xlsxFileNames = cell(Nfiles,1);

for i = 1:Nfiles
    xlsxFileNames{i} = xlsxFiles(i).name;
end

xlsxFileNamesSorted = sort_nat(xlsxFileNames);

for i = 1:Nfiles
    
    sheetName = strrep(xlsxFileNamesSorted{i}, '.xlsx', '');
    
    ImportTable = readtable([xlsxFiles(i).folder,'\',xlsxFileNamesSorted{i}]);
    
    writetable(ImportTable,'SourceData.xlsx','WriteMode','Append','Sheet',sheetName)

end