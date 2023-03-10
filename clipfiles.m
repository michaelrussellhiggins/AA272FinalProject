%% Write Clipped files
close all
clear
clc
%% Load sheet of clip times and indices

clips = readtable('routeclipping.xlsx');
%%

gpsposfiles = dir('Data/GPSPositionDirect2/*.csv');
gpsfilenames = {gpsposfiles.name};

imuposfiles = dir('Data/IMUReadings/*.xlsx');
imufilenames = {imuposfiles.name};

for i = 1:72
    gpsfile = gpsfilenames{i};
    gpsdata = readtable(['Data/GPSPositionDirect/' gpsfile]);
    gpsdataclip = gpsdata(clips.GPSStartIDX(i):clips.GPSStopIDX(i),:);

    imufile = imufilenames{i};
    acceldata = readtable(['Data/IMUReadings/' imufile],'Sheet','accel');
    acceldataclip = acceldata(clips.AccStartIDX(i):clips.AccStopIDX(i),:);
    orientdata = readtable(['Data/IMUReadings/' imufile],'Sheet','orientation');
    orientdataclip = orientdata(clips.OrientStartIDX(i):clips.OrientStopIDX(i),:);

    writeGfile = ['Data/GPSPositionClip2/',gpsfile];
    writetable(gpsdataclip,writeGfile);

    writeIfile = ['Data/IMUReadingsClip/',imufile];
    writetable(acceldataclip,writeIfile,'Sheet','accel');
    writetable(orientdataclip,writeIfile,'Sheet','orientation');

end