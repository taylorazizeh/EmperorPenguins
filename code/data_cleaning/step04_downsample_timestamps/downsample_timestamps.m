% Title: Downsample timestamps
% Author: Taylor Azizeh
% Date: 5 February 2025
% Description: This script is optional - a troubleshooting script if you
% find that you need to downsample timestamps that were at at 100Hz.

% Define input and output folders
timePath = '/Users/taylorazizeh/Documents/Research/active/Emperor penguins/data/cleaned/timestamps/100Hz';
outputPath = '/Users/taylorazizeh/Documents/Research/active/Emperor penguins/data/cleaned/timestamps';

% List all timestamp files
timeFiles = dir(fullfile(timePath, '*_time.nc'));

% Loop through each timestamp file
for i = 1:length(timeFiles)
    % Get the input file path
    inputFile = fullfile(timePath, timeFiles(i).name);
    
    % Read the timestamp data
    time_data = ncread(inputFile, 'Timestamp'); % Replace with actual variable name in .nc files
    
    % Downsample the data by keeping every other row
    time_data_downsampled = time_data(1:2:end, :); % Ensure it remains a character array
    
    % Extract bird ID and create the output file name
    [~, fileName, ~] = fileparts(timeFiles(i).name); % Get the file name without extension
    outputFile = fullfile(outputPath, [fileName, '_50Hz_time.nc']);
    
    % Write the downsampled data to a new .nc file
    nccreate(outputFile, 'Timestamp', ...
        'Dimensions', {'time', size(time_data_downsampled, 1), 'char_dim', size(time_data_downsampled, 2)}, ...
        'Datatype', 'char', ...
        'DeflateLevel', 5); % Enable compression to reduce file size
    ncwrite(outputFile, 'Timestamp', time_data_downsampled);
    
    % Display progress
    fprintf('Processed and saved downsampled timestamp file: %s\n', outputFile);
end