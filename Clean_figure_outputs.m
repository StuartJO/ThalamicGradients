%function Clean_figure_outputs()

% Define the source folder
sourceFolder = './figure_outputs';

% Define the target folders in the source folder
bitmapFolder = fullfile(sourceFolder, 'Bitmap');
vectorFolder = fullfile(sourceFolder, 'Vector');

% Create the target folders if they don't exist
if ~exist(bitmapFolder, 'dir')
    mkdir(bitmapFolder);
end

if ~exist(vectorFolder, 'dir')
    mkdir(vectorFolder);
end

% Process files in the sourceFolder itself
sourceFileList = dir(fullfile(sourceFolder, '*.*'));
for k = 1:length(sourceFileList)
    if ~sourceFileList(k).isdir
        [~, ~, extension] = fileparts(sourceFileList(k).name);
        
        % Move .png files to Bitmap folder
        if strcmpi(extension, '.png')
            sourceFilePath = fullfile(sourceFolder, sourceFileList(k).name);
            targetPath = fullfile(bitmapFolder, sourceFileList(k).name);
            movefile(sourceFilePath, targetPath);
        % Move .svg files to Vector folder
        elseif strcmpi(extension, '.svg')
            sourceFilePath = fullfile(sourceFolder, sourceFileList(k).name);
            targetPath = fullfile(vectorFolder, sourceFileList(k).name);
            movefile(sourceFilePath, targetPath);
        end
    end
end

% Get a list of all folders and subfolders
folderList = dir(sourceFolder);
folderList = folderList([folderList.isdir]);
folderList = folderList(~ismember({folderList.name}, {'.', '..'}));

% Loop through each folder
for i = 1:length(folderList)
    subfolder = fullfile(sourceFolder, folderList(i).name);
    
    % Check if the subfolder is not named "Bitmap" or "Vector"
    if ~strcmpi(folderList(i).name, 'Bitmap') && ~strcmpi(folderList(i).name, 'Vector')
        % Create the target folders if they don't exist in the subfolder
        bitmapFolder = fullfile(subfolder, 'Bitmap');
        vectorFolder = fullfile(subfolder, 'Vector');
        
        if ~exist(bitmapFolder, 'dir')
            mkdir(bitmapFolder);
        end
        
        if ~exist(vectorFolder, 'dir')
            mkdir(vectorFolder);
        end
        
        % Get a list of all files in the subfolder
        fileList = dir(fullfile(subfolder, '*.*'));
        
        % Loop through each file
        for j = 1:length(fileList)
            % Check if it's a file (not a folder)
            if ~fileList(j).isdir
                % Check the file extension
                [~, ~, extension] = fileparts(fileList(j).name);
                
                % Move .png files to Bitmap folder
                if strcmpi(extension, '.png')
                    sourceFilePath = fullfile(subfolder, fileList(j).name);
                    targetPath = fullfile(bitmapFolder, fileList(j).name);
                    movefile(sourceFilePath, targetPath);
                % Move .svg files to Vector folder
                elseif strcmpi(extension, '.svg')
                    sourceFilePath = fullfile(subfolder, fileList(j).name);
                    targetPath = fullfile(vectorFolder, fileList(j).name);
                    movefile(sourceFilePath, targetPath);
                end
            end
        end
    end
end