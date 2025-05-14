%Created by: Rachel Holser (rholser@ucsc.edu), rewritten from P.Robinson's 2009 script. 
%Last updated: 21-Jun-2022

%Function to prepare csv file, fix errors, and run dive analysis on TDR
%data.

%Requires IKNOS toolbox
%Requires functions: yt_iknos_da_RRH
%                    DA_data_compiler_RRH
%                    

function output=change_format_DA2_pinguino_RRH_TV4_0_1(filename)
    %Step 1: load csv of TDR data
    data=readtable(filename,'HeaderLines',0,'ReadVariableNames',true);
    
    %Test for normal headers.  If Depth header is missing, re-imports file
    %with no headers and assigns them.
    HeaderTest=find(strcmp('Depth', data.Properties.VariableNames)==1);
    if isempty(HeaderTest)
        data.Properties.VariableNames={'Time','Depth','Var3','Var4','Var5',...
            'Var6','Var7','Var8'};
    end

    %Step 2: remove rows with depth NaNs and round to nearest 0.1m
    ind_nan=find(isnan(data.Depth));
    data(ind_nan,:)=[];
    
    data.Depth=round(data.Depth,1);

    %Step 3: date conversion to datenum
    try
        data.JDay=datenum(data.Timestamp);
    end

    [data.Year, data.Month, data.Day, data.Hour, data.Min, data.Sec]...
        =datevec(data.JDay);
    data.Sec=round(data.Sec); 

    %Step 4: calculate sampling rate
    SamplingRate=round(mode(diff(data.Sec)));

    %Step 5: generate variable string for iknos da and write to new .csv
    [data_DA,DAstring]=DA_data_compiler_pinguino_RRH_TV4(data);
    writematrix(data_DA,[strtok(filename,'.') '_DAprep.csv']);
    fid=fopen([strtok(filename,'.') '_DAString.txt'],'wt');
    fprintf(fid,DAstring);
    fclose(fid);

    %Step 6: Depth resolution detection
    DepthRes=unique(abs(diff(data.Depth)));
    if DepthRes(1,1)>0
        DepthRes=DepthRes(1,1);
    else
        DepthRes=DepthRes(2,1);
    end

    %Step 7: Run IKNOS DA   
    yt_iknos_da_RRH([strtok(filename,'.') '_DAprep.csv'],DAstring,...
        20/SamplingRate,3/DepthRes,20,'wantfile_yes','ZocWindow',0.5,'ZocMinMax',[-10,800]);
    
end
