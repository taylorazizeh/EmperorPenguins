%Compiles data for input into IKNOS DA and generates corresponding string
%of variables

%updated 26-Feb-2021

function [output, outputstring]=DA_data_compiler_pinguino_RRH_TV4(data)

m=size(data,2);
data_names=data.Properties.VariableNames;

output=[data.Year, data.Month, data.Day, data.Hour, data.Min, data.Sec];
outputstring=['Year Month Day Hour Minute Second'];
for i=1:m
    if strcmp({'Depth'},data_names{1,i}(1:end))==1
        output=[output,data.Depth];
        outputstring=[outputstring ' depth'];
    elseif strcmp({'ExternalTemp'},data_names{1,i}(1:end))==1
        output=[output,data.ExternalTemp];
        outputstring=[outputstring ' etemp'];
    elseif strcmp({'ExternalTemperature'},data_names{1,i}(1:end))==1
        output=[output,data.ExternalTemperature];
        outputstring=[outputstring ' etemp'];
    elseif strcmp({'Temp___C_'},data_names{1,i}(1:end))==1
        output=[output,data.Temp___C_];
        outputstring=[outputstring ' etemp'];
    elseif strcmp({'InternalTemperature'},data_names{1,i}(1:end))==1
        output=[output,data.InternalTemperature];
        outputstring=[outputstring ' itemp'];
    elseif strcmp({'LightLevel'},data_names{1,i}(1:end))==1
        output=[output,data.LightLevel];
        outputstring=[outputstring ' alight'];
    end
end
