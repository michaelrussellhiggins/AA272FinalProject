clc
clear all
close all

% Make sure AA272FinalProjectFolder is open in MATLAB

% Gets list of data files

% files = dir('Data/Sorted/*.txt');
% 
% filenames = {files.name};

% Formats data file names properly

% for i = 1:length(filenames)
%     filenames(i) = erase(filenames(i), '.txt');
% end

% Create seperate GPS and IMU data files from each data file

% for i = 1:numel(filenames)
%     CreateDataFiles('Data/Sorted/',filenames{i}, 'Data/Parsed/');
% end

% Filter non-GPS measurements

% gnssfiles = dir('Data/Parsed/*.csv');
% 
% gnssfilenames = {gnssfiles.name};
% 
% for i = 1:length(gnssfilenames)
%     GetGPSMeasurements(gnssfilenames{i});
% end

gpsfiles = dir('Data/GPSPositionClip/*.csv');

gpsfilenames = {gpsfiles.name};

% Get pseudoranges from GPS data

% for i = 1:length(gpsfilenames)
%     GetPseudoranges(gpsfilenames{i});
% end

% Get GPS positions from ephemeris data

% for i = 1:length(gpsfilenames)
%     GetSatellitePosition(gpsfilenames{i});
% end

% Get user position from GPS data

% for i = 1:length(gpsfilenames)
%     GetUserPosition(gpsfilenames{i});
% end

% Convert user ECEF data to local coordinates



% Extract IMU data from each .xlsx file spreadsheet

imufiles = dir('Data/IMUReadingsClip/*.xlsx');

imufilenames = {imufiles.name};

% for i = 1:length(imufilenames)
%     GetCorrectedIMUData(imufilenames{i});
% end

% Apply Kalman filter

for i = 1:length(imufilenames)
    FuseGPSIMU4(imufilenames{i}(1:end-5))
end

