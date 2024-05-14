function number = getseed(fileName,index)
% get the seed
%% Input
%  fileName: name of file
%  index: the index that you want to extract
%% Output
%  number: gives you the seed number from the txt file named fileName
%%%%%%%%%%%%%%
% Open the file for reading
fileID = fopen(fileName, 'r');

% Check if the file was opened successfully
if fileID == -1
    error('Unable to open file for reading.');
else
    % Read the first entry from the file
    % Skip the first four entries
    for i = 1:index-1
        fscanf(fileID, '%d', 1);
    end
    number = fscanf(fileID, '%d', 1);
    
    % Close the file
    fclose(fileID);
    
end
end