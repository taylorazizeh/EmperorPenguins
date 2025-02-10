% Title: csv to netCDF file conversion
% Author: Taylor Azizeh
% Date: 5 February 2025
% Description: This script extracts variables of your choice from a .csv
% file and exports them as seperate netCDF files.


% THIS IS A LOOP

cd '/Users/taylorazizeh/Library/CloudStorage/GoogleDrive-taylor.azizeh@sjsu.edu/.shortcut-targets-by-id/1qVmvBCQU8aplOO73CUX_mqRay94h7cgI/Penguin_Shared/Data/Crozier/2019/Axytrek/Decimated'


% List of variable names in the CSV and their corresponding names for NetCDF files
varsCsv = {'AccelX', 'Depth', 'JerkX'};
varsNc = varsCsv; % Using the same names for NetCDF files

% List of file identifiers, adjust according to your specific files
fileIds = {'309f'}; % Add the rest of the identifiers

% Loop through each file identifier
for fileId = fileIds
    % Construct the CSV file name
    csvFileName = sprintf('22EP_%s_filtered_acceleration.csv', fileId{1});
    
    % Load the data from the CSV file
    data = readtable(csvFileName);
    
    % Loop through each variable
    for v = 1:length(varsCsv)
        % Construct the NetCDF file name
        ncFileName = sprintf('19EP_%s_50Hz_%s.nc', fileId{1}, varsNc{v});
        
        % Create a NetCDF file for the current variable
        ncid = netcdf.create(ncFileName, 'CLOBBER');
        
        % Define dimensions
        dimid_time = netcdf.defDim(ncid, 'Time', size(data, 1));
        
        % Define variable
        varid = netcdf.defVar(ncid, varsNc{v}, 'double', dimid_time);
        
        % Optional: Define variable attributes
        % You can add appropriate units or other attributes here
        % For example:
        if strcmp(varsNc{v}, 'Depth')
            netcdf.putAtt(ncid, varid, 'units', 'm');
        elseif strcmp(varsNc{v}, 'AccelX')
            netcdf.putAtt(ncid, varid, 'units', 'g');
        elseif strcmp(varsNc{v}, 'JerkX')
            netcdf.putAtt(ncid, varid, 'units', 'g/s');
        end
        
        % End definition mode
        netcdf.endDef(ncid);
        
        % Write data for the current variable
        netcdf.putVar(ncid, varid, table2array(data(:, varsCsv{v})));
        
        % Close the NetCDF file
        netcdf.close(ncid);
    end
end
