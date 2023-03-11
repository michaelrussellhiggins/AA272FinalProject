%ProcessGnssMeasScript.m, script to read GnssLogger output, compute and plot:
% pseudoranges, C/No, and weighted least squares PVT solution
%
close all
clear
clc
%%
gpsfiles = dir('../Data/Sorted/*.txt');

gpsfilenames = {gpsfiles.name};

for i = 1:length(gpsfilenames)
% you can run the data in pseudoranges log files provided for you: 
prFileName = gpsfilenames{i}; %with duty cycling, no carrier phase
% prFileName = 'pseudoranges_log_2016_08_22_14_45_50.txt'; %no duty cycling, with carrier phase
% as follows
% 1) copy everything from GitHub google/gps-measurement-tools/ to 
%    a local directory on your machine
% 2) change 'dirName = ...' to match the local directory you are using:
dirName = '../Data/Sorted/';
% 3) run ProcessGnssMeasScript.m script file 
param.llaTrueDegDegM = [];
loc = 0;
loc_char = prFileName(9:11);
if (strcmp(loc_char,'Out') == 1) || (strcmp(loc_char,'Dra') == 1)
    loc = 1;
elseif (strcmp(loc_char,'Cur') == 1) || (strcmp(loc_char,'Fly') == 1)
    loc = 2;
elseif (strcmp(loc_char,'Cor') == 1) || (strcmp(loc_char,'Pos') == 1)
    loc = 3;
end
%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

%% data

%To add your own data:
% save data from GnssLogger App, and edit dirName and prFileName appropriately
%dirName = 'put the full path for your directory here';
%prFileName = 'put the pseuoranges log file name here';

%% parameters
%param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
%param.llaTrueDegDegM = [37.422578, -122.081678, -28];%Charleston Park Test Site

%% Set the data filter and Read log file
dataFilter = SetDataFilter;
[gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
if isempty(gnssRaw), return, end

%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
utcTime = Gps2Utc([],fctSeconds);
allGpsEph = GetNasaHourlyEphemeris(utcTime,'../Ephem');
if isempty(allGpsEph), return, end

%% process raw measurements, compute pseudoranges:
[gnssMeas] = ProcessGnssMeas(gnssRaw);

%% plot pseudoranges and pseudorange rates
% h1 = figure;
% [colors] = PlotPseudoranges(gnssMeas,prFileName);
% h2 = figure;
% PlotPseudorangeRates(gnssMeas,prFileName,colors);
% h3 = figure;
% PlotCno(gnssMeas,prFileName,colors);

%% compute WLS position and velocity
gpsPvt = GpsWlsPvt2(gnssMeas,allGpsEph);

%% write csv files
utc = gnssMeas.utc(:,1);
pos = gpsPvt.allECEF;
lla = gpsPvt.allLlaDegDegM;
enu = zeros(length(pos),3);

switch loc
    case 1
        lat_station = 37.430779;
        lon_station = -122.155025;
    case 2
        lat_station = 37.431276;
        lon_station = -122.155041;

    case 3
        lat_station = 37.431101;
        lon_station = -122.155149;
end

for j = 1:length(pos)
    r_ecef = pos(j,:);
    [E, N, U] = geodetic2enu(lla(j,1), lla(j,2), lla(j,3), lat_station, lon_station, 16.5,wgs84Ellipsoid);
    enu(j,:) = [E N U];
    %enu(j,:) = (ecefTOenu(r_ecef',loc))';
end

user_pos_mat = [utc pos lla enu];
user_pos_tab = array2table(user_pos_mat,...
    'VariableNames',{'utc','X','Y','Z','lat','lon','h','E','N','U'});

writefile = ['../Data/GPSPositionDirect2/',prFileName(1:end-4),'.csv'];
writetable(user_pos_tab,writefile)

end
%% plot Pvt results
% h4 = figure;
% ts = 'Raw Pseudoranges, Weighted Least Squares solution';
% PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
% h5 = figure;
% PlotPvtStates(gpsPvt,prFileName);

%% Plot Accumulated Delta Range 
% if any(any(isfinite(gnssMeas.AdrM) & gnssMeas.AdrM~=0))
%     [gnssMeas]= ProcessAdr(gnssMeas);
%     h6 = figure;
%     PlotAdr(gnssMeas,prFileName,colors);
%     [adrResid]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
%     h7 = figure;
%     PlotAdrResids(adrResid,gnssMeas,prFileName,colors);
% end
%% end of ProcessGnssMeasScript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016 Google Inc.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
