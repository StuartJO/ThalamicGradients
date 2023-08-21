function deleteSvgLines(filename)
    % Read the .svg file and store its content in a cell array
    fid = fopen(filename, 'r');
    svgContent = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    svgContent = svgContent{1};

    % Process the content to find lines starting with "><circle"
    linesToDelete = find(strncmp(svgContent, '><circle', 8));

    % Remove the lines matching the criteria and the lines before and after them
    linesToDelete = unique(sort([linesToDelete - 1; linesToDelete; linesToDelete + 1; linesToDelete + 2; linesToDelete + 3; linesToDelete + 4]));
    svgContent(linesToDelete) = [];

    % Write the modified content back to the .svg file
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\n', svgContent{:});
    fclose(fid);
end