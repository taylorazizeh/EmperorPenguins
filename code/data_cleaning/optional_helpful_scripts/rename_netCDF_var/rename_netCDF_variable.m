% Title: Rename netCDF variable
% Author: Taylor Azizeh
% Date: 5 February 2025
% Description: This script is optional - a troubleshooting script if you
% need to rename a variable in a netCDF file.

% Define the file path for the problematic bird
filePath = '/Users/taylorazizeh/Documents/Research/active/Emperor penguins/data/cleaned/acceleration/19EP_304f_50Hz_Z.nc'; % Update this to the actual file path

% Temporary file to save the modified NetCDF
tempFilePath = [filePath, '.tmp'];

% Open the original NetCDF file in read mode
ncid_in = netcdf.open(filePath, 'NC_NOWRITE');

% Create a new NetCDF file for output
ncid_out = netcdf.create(tempFilePath, 'CLOBBER');

% Get information about the original NetCDF file
info = ncinfo(filePath);

% Copy dimensions
for i = 1:length(info.Dimensions)
    dimName = info.Dimensions(i).Name;
    dimLength = info.Dimensions(i).Length;
    netcdf.defDim(ncid_out, dimName, dimLength);
end

% Copy variables, renaming the target variable
for i = 1:length(info.Variables)
    varName = info.Variables(i).Name;
    varType = info.Variables(i).Datatype;
    varDims = {info.Variables(i).Dimensions.Name};
    varDimsID = cellfun(@(x) netcdf.inqDimID(ncid_out, x), varDims);

    % Check if this is the variable to rename
    if strcmp(varName, 'Z_acceleration')
        varName = 'Z'; % Rename the variable
    end

    % Define the variable in the new file
    netcdf.defVar(ncid_out, varName, varType, varDimsID);

    % Copy variable attributes
    for j = 1:length(info.Variables(i).Attributes)
        attrName = info.Variables(i).Attributes(j).Name;
        attrValue = info.Variables(i).Attributes(j).Value;
        netcdf.putAtt(ncid_out, netcdf.inqVarID(ncid_out, varName), attrName, attrValue);
    end
end

% End define mode for the new file
netcdf.endDef(ncid_out);

% Copy data for each variable
for i = 1:length(info.Variables)
    varName = info.Variables(i).Name;
    data = ncread(filePath, varName);

    % Write data to the new file, renaming the variable if needed
    if strcmp(varName, 'Z_acceleration')
        netcdf.putVar(ncid_out, netcdf.inqVarID(ncid_out, 'Z'), data);
    else
        netcdf.putVar(ncid_out, netcdf.inqVarID(ncid_out, varName), data);
    end
end

% Close the files
netcdf.close(ncid_in);
netcdf.close(ncid_out);

% Replace the original file with the modified one
movefile(tempFilePath, filePath);

disp('Variable renamed successfully.');