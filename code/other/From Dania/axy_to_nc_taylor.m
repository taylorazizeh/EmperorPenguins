clear all
close all

%% add tagtools to path
addpath(genpath('C:\Users\taylo\OneDrive\Documents\MLML\VEL\Thesis\Data Processing\TagTools'))

%% Specify deployment ID and path to csv file for that depid
depid = 'af19_304f';
fname = 'C:\Users\taylo\OneDrive\Documents\MLML\VEL\Thesis\Data Processing\Data\Raw\2019\100Hz\19EP_304f.csv';

%% Convert the csv file to an nc file
ncfile = read_axytrek_taylor2(fname,depid,'af','bm');

%% Load the new nc file and look at its metadata
% ncfile = 'af19_304f_raw';
load_nc(ncfile)
info

%% Add additional deployment and project information to the info metadata structure
info.dephist_device_tzone = '0';
info.dephist_deploy_locality  = 'Cape Crozier, Antarctica';
info.dephist_deploy_location_lat = -77.26698;
info.dephist_deploy_location_lon = 169.15304;
% info.dephist_deploy_datetime_start = '16-11-2018 01:00:00';
info.dephist_deploy_method = 'tape, cable ties' ;
info.project_name = 'Foraging ecology of late chick-rearing emperor penguins';
info.project_datetime_start = '';
info.project_datetime_end = '';
info = orderfields(info) ;

%% Plot dive profile
figure, plott(P)

%% Plot temperature
figure, plott(T)

%% It seems that the transients in the dive profile correspond in time to temperature errors (extremely low temperatures), so we can perhaps use those to clean up the pressure data
find(T.data<-100)
find(T.data<-100)/T.sampling_rate
P.sampling_rate*(find(T.data<-100)/T.sampling_rate)
P.data(P.sampling_rate*(find(T.data<-100)/T.sampling_rate))
max(P.data)

diff(P.sampling_rate*(find(T.data<-100)/T.sampling_rate)) % single values, so could just assume the previous value, or interpolate between 2 neighbouring values
P.data(P.sampling_rate*(find(T.data<-100)/T.sampling_rate)+(-1:2:1))
mean(P.data(P.sampling_rate*(find(T.data<-100)/T.sampling_rate)+(-1:2:1)),2)

%% Replace the pressure errors with interpolated values
P.data(P.sampling_rate*(find(T.data<-100)/T.sampling_rate)) = mean(P.data(P.sampling_rate*(find(T.data<-100)/T.sampling_rate)+(-1:2:1)),2);
figure, plott(P)

%% There seem to also be some incredibly high acceleration values, some transient, some extending over longer intervals (e.g. on the third axis):
% figure, plott(A)
ix = find(abs(A.data(:,3))>1000*prctile(abs(A.data(:,3)),99.99))
A.data(ix,3)
figure, plot(diff(ix))

%% But those do not correspond to the temperature transients/errors
 A.data(round(A.sampling_rate*(find(T.data<-100)/T.sampling_rate))+(-10:10),:)

%% Rather than interpolating, for now, I'll replace the values with NaNs. This should be kept in mind in future analysis and if causing issues, resolved differently!
for i = 1:3
    ix = find(A.data(:,i)>1000*prctile(A.data(:,i),99.99));
    A.data(ix,i) = NaN;
end

figure, plott(P, A) % the x axes of the subplots are linked, so when zooming on one, the other follows

%% Plot track data
figure, plot(POS.data(:,3),POS.data(:,2),'r.-'), grid on % column names stored in POS.column_name

%% Plot track data with google maps 
figure, plot(POS.data(:,3),POS.data(:,2),'g.-') 
plot_google_map('MapType','hybrid')

%% Compute jerk and generate a sensor structure for it
J = njerk(A);
J = sens_struct(J, A.sampling_rate, depid, 'jerk')

%% Plot dive profile along with jerk
figure, plott(P, J)

%% Fix pressure offset at the surface - you can play with different time intervals and see if the data look better, and if so, if you want to store those values instead below (here, I don't do that)
[pp, poffs] = fix_offset_pressure(P, 15*60, 8*3600);
figure, plott(P, pp) % plot for comparison
% ylim([-5,5])

%% Save an updated nc file with (some of) the processed data
upd_ncfile = [depid, '_sens']; % or whatever you want to call it
save_nc(upd_ncfile, info, A, P, POS, T, J) ;
