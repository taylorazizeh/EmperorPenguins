%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rachel Holser (rholser@ucsc.edu), modified from P.Robinson's code
%Last Updated: 21-Apr-2022
%
%   Compiles tracking, diving, foraging success, and metadata into tables 
%   and structures in a single .mat file.  
%   Uses foieGras processed tracking data to generate Track_Best and
%   find locations for each dive (DiveLoc_Best). Creates _TV4_alpha.mat
%   file.
%
%
%Required functions:  yt_interpol_linear_2
%
%Preparation:
%      Process and assemble all tracking and diving data (Steps 1-3)
%      Create metadata.mat (start/stop, foraging success, tag metadata)

%IMPORTANT: Data filenames MUST use standard endings to be discovered
%       and imported correctly in this script: *_foieGras_crw.csv, 
%       *_argos_raw.csv, *_argos_crawl.csv, *_FastGPS.csv, *_llgeo_raw.csv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%load metadata information for all elephant seal deployments
load MetaData.mat

%Create row of NaNs in TagMetaDataAll structure
TagMetaDataAll=[array2table(nan(1,20),'variablenames',...
    TagMetaDataAll.Properties.VariableNames);TagMetaDataAll];

%find all dive files to be processed
files=dir('*_full.csv');

%find indices of first and last deployments to be processed within MetaData
%structures.  All deployments inbetween the first and last will be
%"compiled", regardless of whether data are present
first=find(MetaDataAll.TOPPID==str2double(strtok(files(1).name,'_')));
last=find(MetaDataAll.TOPPID==str2double(strtok(files(end).name,'_')));
clear files

%Loop through deployments from first to last as defined aboce
for i=first:last 
 %% Step 1: clear variables and create empty variables for new file   
        
    %Display current loop's deployment ID
    disp([MetaDataAll.FieldID{i} '_' num2str(MetaDataAll.TOPPID(i))])
    
    %Create all empty variables to be populated
    TOPPID=MetaDataAll.TOPPID(i);
    DiveStat=array2table([]);
    Track_foieGras=array2table([]);
    DiveLoc_Best=array2table([]);
    MetaData=[];
    
    DepartDate=datenum(MetaDataAll.DepartDate(i));
    ArriveDate=datenum(MetaDataAll.ArriveDate(i));
        
    %% Step 2: MetaData
    
    %%%Populate MetaData%%%
    %Copy from MetaDataAll table and convert to structure
    MetaData=table2struct(MetaDataAll(i,:));
    %Create both datenum and datestring for departure and arrival
    MetaData.DepartDateTime=MetaData.DepartDate;
    MetaData.DepartDate=datenum(MetaData.DepartDate);
    MetaData.ArriveDateTime=MetaData.ArriveDate;
    MetaData.ArriveDate=datenum(MetaData.ArriveDate);
    
    %%%Populate MetaData.TagDeployInfo%%%
    %Find row in TagMetaDataAll that matches ToppID and copy over data
    row=find(TagMetaDataAll.TOPPID(:)==TOPPID);
    %if no matching row is found, populate MetaData.TagDeployInfo with row of NaNs
    if isempty(row)
        MetaData.TagDeployInfo=table2struct(TagMetaDataAll(1,:));       
    else
        MetaData.TagDeployInfo=table2struct(TagMetaDataAll(row,:));
    end
    clear row
    
    %%%Create list of TOPP IDs for other deployments on same animal
    rows=strcmp(MetaData.FieldID,MetaDataAll.FieldID);
    MetaData.AllDeployments=MetaDataAll.TOPPID(rows,1);
    clear rows
    
    %% Step 3: Check for dive of data and read if present, truncate to startstop
    %times.  
    
    %find DiveStat file for current deployment
    divefile=dir([num2str(TOPPID) '*DiveStat.csv']);
    
    if size(divefile,1)==0
        %No dive data, move to Step 4, tracking data
        MetaData.TDRused.DiveType=[];
        MetaData.TDRused.DiveID=[];
    elseif size(divefile,1)==1
        %Single TDR available, no prioritization needed
        DiveStat=readtable(divefile(end).name);
    end
    
    %truncate DiveStat to depart/arrive datetime
    if ~isnan(ArriveDate)==1
        DiveStat=DiveStat(DiveStat.JulDate>=DepartDate &...
            DiveStat.JulDate<=ArriveDate,:);
    else
        DiveStat=DiveStat(DiveStat.JulDate>=DepartDate,:);
    end
     
%% Step 4: check for tracking data and read what is present into tables,
%label variables, truncate to depart/arrive dates

%%%%%IMPORTANT: Data filenames MUST use standard endings to be discovered
%%%%%and imported correctly in this script: *_foieGras_crw.csv, 

    %Check for presence of data type. 
    datafile=dir([num2str(TOPPID) '*_foieGras_crw.csv']);
    if size(datafile,1)>0
        %If file is present, load using readtable
        Track_foieGras=readtable(datafile(end).name);
        %Remove leading variable (numeric index) from foieGras files, if
        %present
        try
            Track_foieGras=removevars(Track_foieGras,{'Var1'});
        end
        %Rename columns for consistency
        Track_foieGras.Properties.VariableNames = {'TOPPID','DateTime','Lon',...
            'Lat','x','y','x_se_km','y_se_km','u','v','u_se_km','v_se_km'};
        %Create matlab JulDate
        Track_foieGras.JulDate=datenum(Track_foieGras.DateTime);
    end

%% Step 6: Create DiveLoc_Best
    
    %Is there diving data?  If not, do nothing.
    if size(divefile,1)==0
        
    elseif size(divefiles,1)==1
        try
            %Linear interpolatation of Track_Best to match dive datestamps
            %interpolates the track time/lat/lon to the time provided in
            %the divestat data - using that latter time as a "timestepping
            %point"
            DiveLatLon = yt_interpol_linear_2(table2array(Track_foieGras(:,{'JulDate','Lat',...
                'Lon'})),DiveStat.JulDate(:));
            DiveLoc_Best=[DiveLatLon(:,2:3) DiveLatLon(:,1)];
            %Convert to table and rename variables
            DiveLoc_Best=array2table(DiveLoc_Best);
            DiveLoc_Best.Properties.VariableNames={'Lat','Lon','JulDate'};
            %Create Lon360 variable for flexible mapping choices
            DiveLoc_Best.Lon360=DiveLoc_Best.Lon+360;
            ind=find(DiveLoc_Best.Lon>0);
            DiveLoc_Best.Lon360(ind)=DiveLoc_Best.Lon360(ind)-360;
            clear DiveLatLon
        end
    end
    
%% Step 7: Foraging Success and Additional MetaData
    
    %%%check if complete TDR record exists%%%
    if isempty(DiveStat)
        MetaData.Group.CompleteTDR=0;
    else
        %Pct days w/dive data is >95% of trip length?
        test1=size(unique(floor(DiveStat.JulDate)),1)>...
            ((datenum(MetaData.ArriveDate)-datenum(MetaData.DepartDate))*0.95);
        %First dive within 2 days of departure date from starstop?
        test2=(DiveStat.JulDate(1)-datenum(MetaData.DepartDate))<2;
        %Last dive within 2 days of arrival date from starstop?
        test3=(datenum(MetaData.ArriveDate)-DiveStat.JulDate(end))<2;
        %All 3 parameters must be met
        if  (test1+test2+test3)==3 
            MetaData.Group.CompleteTDR=1;
        else
            MetaData.Group.CompleteTDR=0;
        end
    end
    clear test1 test2 test3

    %%%check if complete track exists%%
    if isempty(Track_Best)
        MetaData.Group.CompleteTrack=0;
    else
        %Pct days w/track data is >75% of trip length?
        test1=size(unique(floor(Track_Best.JulDate)),1)>...
            ((datenum(MetaData.ArriveDate)-datenum(MetaData.DepartDate))*0.75);
        %First location within 2 days of departure date from starstop?
        test2=(Track_Best.JulDate(1)-datenum(MetaData.DepartDate))<2;
        %Last location within 2 days of arrival date from starstop?
        test3=(datenum(MetaData.ArriveDate)-Track_Best.JulDate(end))<2;
        %Arrival location is known?
        test4=~isempty(MetaData.ArriveLoc);
        %All 4 parameters must be met
        if (test1+test2+test3+test4)==4 
            MetaData.Group.CompleteTrack=1;
        else
            MetaData.Group.CompleteTrack=0;
        end
    end
    clear test1 test2 test3 test4
       
    %%%%%%%%%%%%%%%%%%%%%%%%   
    %Create new mat file
    save([num2str(TOPPID) '_' MetaData.FieldID '_TV4_alpha.mat'],...
        'TOPPID',...
        'DiveStat',...
        'Track_foieGras',...
        'Track_Best',...
        'DiveLoc_Best',...
        'MetaData')
end