clc
clear all
close all

% Make sure AA272FinalProjectFolder is open in MATLAB

% Gets list of data files

files = dir('Data/Sorted/*.txt');

filenames = {files.name};

% Formats data file names properly

for i = 1:length(filenames)
    filenames(i) = erase(filenames(i), '.txt');
end

% Create seperate GPS and IMU data files from each data file

for i = 1:numel(filenames)
    CreateDataFiles('Data/Sorted/',filenames{i}, 'Data/Parsed/');
end

% Extract GPS data from each .csv file



% Extract IMU data from each .xlsx file spreadsheet

